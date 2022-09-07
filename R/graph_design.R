DiagrammeR::grViz('digraph {

  graph [fontsize=10 fontname="Verdana" compound = true rankdir = LR];
  node [shape=record fontsize=10 fontname="Verdana" style = filled];

  data1 [label = "Dataset Y", shape = folder, fillcolor = Beige];

  data2 [label = "Dataset X\n(optional)", shape = folder, fillcolor = Beige];

  data3 [label = "Dictionary", shape = folder, fillcolor = Beige];

  subgraph cluster_0 {
    node [style=filled];
    "Filter and\ncheck input" "Return Web"

    subgraph cluster_mode {
    node [style=filled];
    "Interaction\nmode" "Membership\nmode"
    label = "Generate binary\nadjacency matrix";
    color=blue;
    }

    label = "Create a web object with\nweaveWebFromTables()";
    color=blue;
  }

  data4 [label = "anansiWeb\nobject", shape = folder, fillcolor = Beige];

  subgraph cluster_1 {
    node [style=filled];
    "Check call" "Pairwise\nassociations" "Differential\nassociations"
    subgraph cluster_diff {
    node [style=filled];
    label = "Differential association\ntesting options";
    color=blue;

    subgraph cluster_df {
    node [style=filled];
    label = "Model types";
    color=blue
    "linear model" "linear mixed\neffects model"} ->
    subgraph cluster_dif {
    node [style=filled];
    label = "Model structre";
    color=blue;
    "Classical\nY ~ X" "Log-ratio\nlog(Y/X)"}

    }
    subgraph cluster_ass {
    node [style=filled];
    label = "Pairwise association\ntesting options";
    color=blue;
    "Between dataset\n(correlation)" "Within dataset\n(propr)"}
    label = "Run analysis in the main function\nanansi()";
    color=blue;
"Collect results" "Account for\nFDR" "Return Yarn"  }


  data5 [label = "anansiYarn\nobject", shape = folder, fillcolor = Beige];

  subgraph cluster_2 {
    node [style=filled];
    "spinToWide()" "spinToLong()" "spinToPlots()";
    label = "Parse output to\nhandy formats";
    color=blue;
  }
  {data1 data2, data3} -> "Filter and\ncheck input"

  "Plot all\nassociations"[shape = oval, fillcolor = lightblue]
  "Plot individual\nassociations"[shape = oval, fillcolor = lightblue]
  "Export as table"[shape = oval, fillcolor = lightblue]


  // Edges between nodes render fine
  "Filter and\ncheck input" -> {"Membership\nmode" "Interaction\nmode"} -> "Return Web";
  "Check call" -> {"Pairwise\nassociations" "Differential\nassociations"}

  "Differential\nassociations" -> {"linear model" "linear mixed\neffects model"}
  {"Classical\nY ~ X" "Log-ratio\nlog(Y/X)"} -> "Collect results"
  "Pairwise\nassociations"      -> {"Between dataset\n(correlation)" "Within dataset\n(propr)"} -> "Collect results"

  "Collect results" -> "Account for\nFDR" -> "Return Yarn" -> "data5"

  "spinToLong()"  -> "Plot all\nassociations"
  "spinToPlots()" -> "Plot individual\nassociations"
  "spinToWide()"  -> "Export as table"

  // Edges that directly connect one cluster to another
  "data4" -> "Check call" [lhead=cluster_1];
  "Return Web" -> "data4" [ltail=cluster_0];

 {"data5"} -> {"spinToLong()"} [ltail=cluster_1 lhead=cluster_2];

}')

