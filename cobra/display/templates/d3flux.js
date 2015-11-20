var model = {{ modeljson }};

require.config({
  paths: {
    d3: "http://d3js.org/d3.v3.min",
    math: "https://cdnjs.cloudflare.com/ajax/libs/mathjs/2.4.0/math.min",
    FileSaver: "https://cdnjs.cloudflare.com/ajax/libs/FileSaver.js/2014-11-29/FileSaver.min",
    d3tip: "https://cdnjs.cloudflare.com/ajax/libs/d3-tip/0.6.7/d3-tip"
  }
});

require(["d3", "math", "FileSaver", "d3tip"], function (d3, math, FileSaver, d3tip) {

  function main(model) {

    var width = 1024,
    height = 768;

    var color = d3.scale.category20();

    var rxncolor = d3.scale.ordinal()
      .range(["#bbb", "#d62728", "#ff7f0e", "#2ca02c", "#9467bd",
              "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"])
      .domain([undefined, 'ko', 1, 2, 3, 4, 5, 6, 7, 8]);

    var force = d3.layout.force()
      .linkDistance(30)
      .charge(-100)
      .chargeDistance(400)
      .gravity(.015)
      .size([width, height]);

    if ({{ no_background }}) {
      var svg = d3.select("#{{ figure_id }}").append("svg")
        .attr("width", width)
        .attr("height", height);
    } else {
      var svg = d3.select("#{{ figure_id }}").select("svg");
    }

    d3.select("#{{ figure_id }}_options .reactionbutton").on("click", function() {
      var $this = $(this);
      $this.toggleClass('btn-danger');
      d3.selectAll(".node.rxn")
        .classed("hidden", function (d, i) {
          return !d3.select(this).classed("hidden");
        });
      if($this.hasClass('btn-danger')){
        $this.text('Hide Reaction Nodes');
      } else {
        $this.text('Show Reaction Nodes');
      }
    });

    d3.select("#{{ figure_id }}_options .svgbutton").on("click", function() { 
      // var svg = d3.select("#{{ figure_id }}").select("svg");
      // var image = Pancake(svg[0][0]);
      // image.download("{{ figure_id }}.png");


      var e = document.createElement('script'); 
      if (window.location.protocol === 'https:') { 
        e.setAttribute('src', 'https://rawgit.com/NYTimes/svg-crowbar/gh-pages/svg-crowbar.js'); 
      } else { 
        e.setAttribute('src', 'http://nytimes.github.com/svg-crowbar/svg-crowbar.js'); 
      } 
      e.setAttribute('class', 'svg-crowbar'); 
      document.body.appendChild(e); 
    });

    function calc_imag_angle(x1, y1, x2, y2) {
      var calcAngleVal = Math.atan2(x1-x2,y1-y2);
      return math.exp(math.multiply(math.i, calcAngleVal));
    }

    function average_angles(rxn, rstoich, nodes) {
      var angles = []
      for (var n in rstoich) {
        var angle = calc_imag_angle(rxn.x, rxn.y, nodes[n].x, nodes[n].y);
        angles.push([math.multiply(angle, rstoich[n])]);
      }
      return math.mean(angles).toPolar().phi;
    }

    function average_dist(rxn, rstoich, nodes) {
      var dist = []
      for (var n in rstoich) {
        var this_dist = math.sqrt(math.square(rxn.x - nodes[n].x) 
            + math.square(rxn.y - nodes[n].y));
        dist.push([this_dist]);
      }
      return math.mean(dist);
    }

    function calculate_path(d, force) {
      var s = d.source,
      t = d.target,
      r = d.rxn,
      cp = {};

      var angle = average_angles(r, d.rstoich, force.nodes());
      var dist = average_dist(r, d.rstoich, force.nodes());

      cp["x"] = r.x - math.multiply(.5*dist, math.sin(angle));
      cp["y"] = r.y - math.multiply(.5*dist, math.cos(angle));

      // return "M" + s.x + "," + s.y + " L " + t.x + "," + t.y;
      // return "M" + s.x + "," + s.y
      //   +" C" + s.x + "," + s.y
      //   + " " + cp.x + "," + cp.y
      //   + " " + r.x + "," + r.y
      //   +" S" + t.x + "," + t.y
      //   + " " + t.x + "," + t.y;
      return "M" + s.x + "," + s.y
        + " Q" + cp.x + "," + cp.y
        + " " + r.x + "," + r.y
        +" T" + t.x + "," + t.y
    }

    var metabolites = jQuery.extend(true, [], model.metabolites),
    reactions = jQuery.extend(true, [], model.reactions),
    nodes = [],
    node_lookup = {},
    links = [],
    mlinks = [],
    bilinks = []
    rxn_stoich = {};

    var mfluxes = [];
    metabolites.forEach(function(metabolite) {
      // Don't add hidden reactions
      if ('notes' in metabolite) {
        if ('map_info' in metabolite.notes) {
          if (metabolite.notes.map_info.hidden) {
            return;
          } else {
            if ({{ hide_unused }} && (Math.abs(metabolite.notes.map_info.flux) < 1E-6)) {
              return;
            }
            mfluxes.push(metabolite.notes.map_info.flux);
          }
        }
      }
      // It's not hidden, add it to the nodes
      nodes.push(metabolite);
    });

    var labelAnchors = [];
    var labelAnchorLinks = [];

    for (var i = 0, len = nodes.length; i < len; i++) {
      var node = nodes[i]

      node_lookup[node.id] = i;
      labelAnchors.push({
        node : node
      });
      labelAnchors.push({
        node : node
      });
      labelAnchorLinks.push({
        source : i * 2,
        target : i * 2 + 1,
        weight : 1
      });
    }

    var fluxes = []
    reactions.forEach(function(reaction) {

      // Don't add hidden reactions
      if ('notes' in reaction) {
        if ('map_info' in reaction.notes) {
          if (reaction.notes.map_info.hidden) {
            return;
          } else if ('flux' in reaction.notes.map_info) {
            if ({{ hide_unused }} && (Math.abs(reaction.notes.map_info.flux) < 1E-6)) {
              return;
            }
            fluxes.push(Math.abs(reaction.notes.map_info.flux));
            if (reaction.notes.map_info.flux < 0) {
              // If the reaction is flowing in reverse, switch products and
              // reactants.
              for (var item in reaction.metabolites) {
                reaction.metabolites[item] *= -1;
              }
            }
          }
        }
      }

      reaction['reactants'] = []
      reaction['products'] = []

      for (var item in reaction.metabolites) {
        if (reaction.metabolites[item] > 0) {
          if (item in node_lookup) {
            // Only add if the node hasn't been hidden
            reaction.products.push(item);
          }
        } else if (item in node_lookup) {
          reaction.reactants.push(item);
        }
      }

      var r_length = reaction.reactants.length,
      p_length = reaction.products.length,
      r_node = {
        "id" : reaction.id,
        "type" : "rxn"
      };

      // Add notes to reaction, if it exists (for map_info)
      if ("notes" in reaction) {
        r_node["notes"] = reaction.notes;
      }

      if (r_length == 0 || p_length == 0) {
        return; // Don't add links on the boundary
      }

      nodes.push(r_node)
        rindex = nodes.length - 1


        if (r_length >= p_length) {
          reaction.reactants.forEach(function (reactant, i) {
            mlinks.push({
              "source" : node_lookup[reactant],
              "target" : node_lookup[reaction.products[i % p_length]],
              "rxn" : rindex
            });
          });
        } else {
          reaction.products.forEach(function (product, i) {
            mlinks.push({
              "source" : node_lookup[reaction.reactants[i % r_length]],
              "target" : node_lookup[product],
              "rxn" : rindex
            });
          });
        }
    });

    mlinks.forEach(function(link) {
      if (link.rxn in rxn_stoich) {
        rxn_stoich[link.rxn][link.source] = 1;
        rxn_stoich[link.rxn][link.target] = -1;
      } else {
        rxn_stoich[link.rxn] = {};
        rxn_stoich[link.rxn][link.source] = 1;
        rxn_stoich[link.rxn][link.target] = -1;
      }
    });

    mlinks.forEach(function(link) {
      var s = nodes[link.source],
      t = nodes[link.target],
      r = nodes[link.rxn];

      links.push({source: s, target: r}, {source: r, target: t});
      bilinks.push({
        "source" : s,
        "target" : t,
        "rxn" : r,
        "rstoich" : rxn_stoich[link.rxn],
      });
    });

    nodes.forEach( function (node) {
      if ("notes" in node) {
        if ("map_info" in node.notes) {
          if (("x" in node.notes.map_info) && ("y" in node.notes.map_info)) {
            node.x = node.notes.map_info.x;
            node.y = node.notes.map_info.y;
            node.fixed = 1;
          }
        }
      }
    });

    // Modify link strength based on flux:
    // link_strength_scale = d3.scale.pow().exponent(1/2)
    link_strength_scale = d3.scale.linear()
      .domain([d3.min(fluxes), d3.max(fluxes)])
      .range([.2, 2]);

    force
      .linkStrength(function (link) {
        try {
          return link.rxn.notes.map_info.flux;
        }
        catch(err) {
          return 1;
        }
      });

    force
      .nodes(nodes)
      .links(links)
      .start();

    svg.append("defs").selectAll("marker")
      .data(reactions)
      .enter()
      .append("marker")
      .attr("id", function (d) { return d.id; })
      .attr("viewBox", "0 0 10 10")
      .attr("refX", 12)
      .attr("refY", 5)
      .attr("markerUnits", "userSpaceOnUse")
      .attr("markerWidth", "8pt")
      .attr("markerHeight", "8pt")
      .attr("orient", "auto")
      .attr("class", "endmarker")
      .append("path")
      .attr("d", "M 0 0 L 10 5 L 0 10 z");

    var link = svg.append('g').selectAll(".link")
      .data(bilinks)
      .enter()
      .append("path")
      .attr("class", function (d) { return "link " + d.rxn.id; })
      .attr("marker-end", function(d) {
        return "url(#" + d.rxn.id + ")"; 
      });


    var node_drag = force.drag()
      .on("dragstart", dragstart);

    function dragstart(d) {
      d3.select(this).classed("fixed", d.fixed = true);
    }

    function releasenode(d) {
      // of course set the node to fixed so the force doesn't include the node in
      // its auto positioning stuff
      d.fixed = false; 
      force.resume();
    }

    // define the nodes
    var node = svg.append("g").selectAll(".node")
      .data(nodes)
      .enter()
      .append("g")
      .on('dblclick', releasenode)
      .call(node_drag);


    node.append("circle")
      .attr("class", function(d) {
        if (d.type == 'rxn') {
          return "node " + d.type + " hidden";
        } else return "node metabolite";})
      .attr("id", function(d) { return d.id; })
      .attr("r", 5)
      .style("fill", function(d) { 
        if (d.type != 'rxn') {
          return color(d.group);
        } else return "";});

    // add the text 
    node.append("text")
      .attr("class", "nodelabel")
      .attr("dx", "12")
      .attr("dy", ".35em")
      .text(function(d) { return d.name; });

    var updateNode = function() {
      this.attr("transform", function(d) {
        return "translate(" + d.x + "," + d.y + ")";
      });
    }

    var updateAnchorLink = function() {
      this.attr("x1", function(d) {
        return d.source.x;
      }).attr("y1", function(d) {
        return d.source.y;
      }).attr("x2", function(d) {
        return d.target.x;
      }).attr("y2", function(d) {
        return d.target.y;
      });
    }

    var updateLink = function() {
      this.attr("d", function(d) {
        return calculate_path(d, force);
      });
    }

    var updateAnchorNodes = function() {
      this.each(function(d, i) {
        if(i % 2 == 0) {
          d.x = d.node.x;
          d.y = d.node.y;
        } else {
          var b = this.childNodes[1].getBBox();

          var diffX = d.x - d.node.x;
          var diffY = d.y - d.node.y;

          var dist = Math.sqrt(diffX * diffX + diffY * diffY);

          var shiftX = b.width * (diffX - dist) / (dist * 2);
          shiftX = Math.max(-b.width, Math.min(0, shiftX));
          var shiftY = 5;
          this.childNodes[1].setAttribute("transform", "translate(" + shiftX + "," + shiftY + ")");
        }
      });
    }

    force.on("tick", function() {
      link.call(updateLink);
      node.call(updateNode);
    });

    // flux_scale = d3.scale.pow().exponent(1/2)
    flux_scale = d3.scale.linear()
      .domain([d3.min(fluxes), d3.max(fluxes)])
      .range([1.5, 8]);

    // metabolite_scale = d3.scale.pow().exponent(1/2)
    metabolite_scale = d3.scale.linear()
      .domain([d3.min(mfluxes), d3.max(mfluxes)])
      .range([4, 8]);

    arrowhead_scale = d3.scale.linear()
      .domain([1.5, 6])
      .range([7, 14]);

    function get_flux_width (rxn) {
      try {
        var flux = flux_scale(Math.abs(rxn.notes.map_info.flux));
        if (!isNaN(flux)) {
          return flux;
        } else{
          return 2.0;
        }
      }
      catch(err) {
        return 2.0; // Default linewidth
      }
    }

    function get_flux_dasharray (d) {
      try {
        if (d.rxn.notes.map_info.group == 'ko') {
          return "5, 5, 1, 5";
        }
        else if (Math.abs(d.rxn.notes.map_info.flux) < 1E-6) {
          return "5,5";
        }
      }
      catch(err) {
        return;
      }
    }

    function get_flux_stroke (rxn) {
      // For now, just keep the same color
      return rxncolor(rxn.notes.map_info.group);
      // return 
    }

    function markerscale (d) {
      return arrowhead_scale(get_flux_width(d)) + "pt";
    }

    function get_node_radius (d) {
      try {
        var nodewidth = metabolite_scale(d.notes.map_info.flux);
        if (!isNaN(nodewidth)) {
          return nodewidth;
        } else {
          return 5;
        }
      }
      catch(err){ return 5; }
    }

    svg.selectAll(".link")
      .attr("stroke-width", function (d) {return get_flux_width(d.rxn);})
      .attr("stroke", function (d) {return get_flux_stroke(d.rxn);})
      .attr("stroke-dasharray", get_flux_dasharray);

    svg.selectAll("marker")
      .attr("markerWidth", markerscale)
      .attr("markerHeight", markerscale)
      .select("path")
      .attr("fill", get_flux_stroke);

    svg.selectAll(".metabolite")
      .attr("r", get_node_radius);

    d3.select("#{{ figure_id }}_options .download")
      .on("click", function () {

        // Add position data to model nodes
        force.nodes().forEach(function (node) {
          if (node.fixed) {
            if (node.type == "rxn") {
              // Reaction object
              obj = $.grep(model.reactions, function(e){ return e.id == node.id; })[0];
              if (!("notes" in obj)) { obj.notes = {}; }
              obj.notes['map_info'] = {
                "x" : node.x,
                "y" : node.y,
              };
            } else {
              // Look in metabolites
              obj = $.grep(model.metabolites, function(e){ return e.id == node.id; })[0];
              if (!("notes" in obj)) { obj.notes = {}; }
              obj.notes['map_info'] = {
                "x" : node.x,
                "y" : node.y,
              };
            }
          }
        });

        var json = JSON.stringify(model);
        var blob = new Blob([json], {type: "application/json"});
        saveAs(blob, "model.json");
      });
  }

  main(model);
});
