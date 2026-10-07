(* ::Package:: *)
(* LieART timings for the cases in benchmark/benchmarks.jl.

   Run from a terminal (LieART must be installed, e.g. via PacletInstall):

       wolframscript -file benchmark/lieart_benchmark.wl > lieart_results.csv

   Each case is run once (after a small warm-up), cut off at 300 s.  Output is
   CSV: task, case, irrep dimension(s), number of results, seconds.  The
   dimensions let us check that LieART's Dynkin-label conventions match Kiwi's
   (they can differ for the exceptional algebras).  Lines starting with "#" are
   diagnostics.
*)

Needs["LieART`"];

$timeLimit = 300;

(* Kiwi (algebra, rank) -> LieART series symbol *)
series = <|"A" -> A, "B" -> B, "C" -> C, "D" -> D, "E" -> E, "F" -> F, "G" -> G|>;
irrep[s_String, labels_List] := Irrep[series[s]] @@ labels;

time[expr_] := Module[{r},
  r = TimeConstrained[AbsoluteTiming[expr], $timeLimit, $Failed];
  If[r === $Failed, {"> " <> ToString[$timeLimit], $Failed}, {NumberForm[First[r], {8, 4}], Last[r]}]
];
SetAttributes[time, HoldAll];

(* Number of distinct irreps in a LieART decomposition (IrrepPlus / IrrepTimes) *)
ncomp[$Failed] := "-";
ncomp[x_] := Length[DeleteDuplicates[Cases[{x}, _Irrep, Infinity]]];

show[s_, l_] := s <> ToString[Length[l]] <> " " <> StringJoin[Riffle[ToString /@ l, ","]];
emit[fields__] := Print[StringRiffle[ToString /@ {fields}, ","]];

characterCases = {
  {"A", {10, 10}}, {"A", {2, 1, 1, 2}}, {"D", {1, 1, 0, 1, 1}}, {"C", {1, 1, 1, 1}},
  {"G", {6, 6}}, {"F", {1, 1, 0, 1}}, {"E", {1, 1, 0, 0, 1, 1}},
  {"E", {0, 0, 1, 0, 0, 0, 1}}, {"E", {1, 0, 0, 0, 0, 0, 1, 0}}
};

tensorCases = {
  {"A", {8, 5}, {6, 7}}, {"A", {1, 1, 1, 1}, {2, 1, 0, 1}},
  {"D", {0, 1, 0, 1, 1}, {1, 0, 1, 0, 0}}, {"G", {3, 3}, {4, 2}},
  {"F", {0, 0, 1, 1}, {1, 0, 0, 1}}, {"E", {1, 1, 0, 0, 0, 1}, {0, 1, 0, 0, 1, 1}},
  {"E", {0, 0, 0, 1, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0}},
  {"E", {0, 0, 0, 0, 0, 0, 1, 1}, {1, 0, 0, 0, 0, 0, 1, 0}}
};

(* warm-up *)
DecomposeProduct[Irrep[A][1], Irrep[A][1]];
WeightSystem[Irrep[A][1, 1]];

emit["task", "case", "dim", "n_results", "seconds"];

(* Characters: full weight system (Kiwi's `character`).  n_results = number of
   weights returned by WeightSystem. *)
Do[
  With[{ir = irrep[c[[1]], c[[2]]]},
    Module[{t = time[WeightSystem[ir]]},
      emit["character", show[c[[1]], c[[2]]], Dim[ir],
        If[t[[2]] === $Failed, "-", Length[t[[2]]]], t[[1]]]]],
  {c, characterCases}];

(* Tensor products *)
Do[
  With[{x = irrep[c[[1]], c[[2]]], y = irrep[c[[1]], c[[3]]]},
    Module[{t = time[DecomposeProduct[x, y]]},
      emit["tensor", show[c[[1]], c[[2]]] <> " x " <> StringJoin[Riffle[ToString /@ c[[3]], ","]],
        ToString[Dim[x]] <> "x" <> ToString[Dim[y]], ncomp[t[[2]]], t[[1]]]]],
  {c, tensorCases}];

(* Plethysms: LieART's API for these is version dependent, so report what exists. *)
Print["# LieART version: ", ToString[Quiet[LieART`$Version]]];
Print["# candidate plethysm functions: ",
  Select[Names["LieART`*"], StringContainsQ[#, "Sym" | "Alt" | "Pleth" | "Tensor", IgnoreCase -> True] &]];
