DiagrammeR::grViz('digraph G {

  graph [fontsize=10 fontname="Verdana" compound=true rankdir = LR];
  node [shape=record fontsize=10 fontname="Verdana" style = filled];

  data1 [label = "TableY", shape = folder, fillcolor = Beige];

  data2 [label = "TableX", shape = folder, fillcolor = Beige];

  data3 [label = "Dictionary", shape = folder, fillcolor = Beige];

  subgraph cluster_0 {
    node [style=filled];
    "Filter" "Create Web";
    label = "Create a web object with\nweaveWebFromTables()";
    color=blue;
  }

  subgraph cluster_1 {
    node [style=filled];
    "Check call" "Pairwise\nassociations" "Differential\nassociations";
    label = "Run analysis in the main function\nanansi()";
    color=blue;
  }

  subgraph cluster_2 {
    node [style=filled];
    "spinToWide()" "spinToLong()" "spinToPlots()";
    label = "Parse output to\nhandy formats";
    color=blue;
  }
  {data1 data2, data3} -> "Filter"

  // Edges between nodes render fine
  "Filter" -> "Create Web";
  "Check call" -> {"Pairwise\nassociations" "Differential\nassociations"}

  "spinToLong()"  -> "Plot all\nassociations"
  "spinToPlots()" -> "Plot individual\nassociations"
  "spinToWide()"  -> "Export as table"

  // Edges that directly connect one cluster to another
  "Create Web" -> "Check call" [ltail=cluster_0 lhead=cluster_1];
 {"Pairwise\nassociations" "Differential\nassociations"} ->
 {"spinToLong()"} [ltail=cluster_1 lhead=cluster_2];
}')

