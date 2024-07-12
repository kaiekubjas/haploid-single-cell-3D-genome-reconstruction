(* ::Package:: *)

(* This file shows some Mathematica code for reading and aggregating real data.

 * Note that Mathematica is needed for running the code.
 * Furthermore, an input text file with lines of the form "#chr_A	pos_A	chr_B	pos_B" with a header is needed.
 *)


(* import the file filename and convert it *)
get=StringSplit[#,"	"]&/@Rest[StringSplit[Import["path/filname","Text"],"\n"]]


(* select those lines with A=B=1 *)
get1=Select[Select[get,First[#]==#[[3]]&],First[#]=="chr1"&]


(* read contacts as edges and beads as vertices 
 * contract 1000000 consecutive vertices, interpreting them as one
 * in the resulting edges, multiple entries are ignored and loops are deleted
 *)
contacts1=Ceiling[(ToExpression[get1[[1;;-1,{2,4}]]]-2999999)/1000000];
contacts2=DeleteDuplicates[contacts1];
contacts3=Select[contacts2,#[[1]]!=#[[2]]&];
Length[contacts3]


(* check whether the resulting graph is connected *)
cgraph=Graph[contacts3];
ConnectedGraphQ[cgraph]


(* show the graph and its backbone (i.e. the path of consecutive vertices) *)
HighlightGraph[cgraph,UndirectedEdge@@#&/@Partition[Range[VertexCount[cgraph]],2,1]]


(* generate the adjacency matrix in a matlab style to be used as a HiC-matrix *)
GenerateMatlabHiC[edges_List]:=StringReplace[StringReplace[ToString[Normal[AdjacencyMatrix[Graph[edges]]]],"}, {"->"\n"],{"{"->"","}"->"",", "->"\t"}]

mat=GenerateMatlabHiC[contacts3];

(* export the file to "folder" with "filename" *)
Export["folder/filename.txt",mat]
