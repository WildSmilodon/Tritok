*#
*# GiD 7.2 script za generiranje 3D mreze
*# spisal Jure Ravnik, Februar 2007
*#
*# Uporaba :
*#
*# - v GiDu izberi problem type = bem3d
*# - narisi geometrijo iz crt
*# - naredi NUBS surface
*# - izberi Data/Conditions
*#          - tu povej katera crta ima katero stevilko stene : od 1 naprej
*#            lahko imata dve skupaj dotikajoci se crti isto cifro
*#
*# - Generiranje mreze
*#   Meshing/Quadratic elements/Quadratic9
*#   Meshing/Structured/Surface : doloci obmocja, nato za vsako crto stevilo elementov
*#   Meshing/generate
*#
*# - Zazeni izpis *.dat datoteke s Calculate/calcualte
*#
# Version
VERS GiD-3D-1.0

*# ------------------------------------------
*# Prestejem robna vozlisca
*# ------------------------------------------
*Set Cond boundary *nodes
*loop nodes *OnlyInCond
*set var NBNODES=LoopVar
*end

# Title
*# ------------------------------------------
*# Dolocim npo, npoc
*# ------------------------------------------
*#if((nnode==9))
*set var NPOB=9
Qua.el.
*#endif

*# ------------------------------------------
*# Dolocim nbelem
*# ------------------------------------------
*Set Cond boundary *elems *CanRepeat
*loop elems *OnlyInCond
*set var NBELEM=LoopVar
*end
*# ------------------------------------------
*# izpis velikosti
*# ------------------------------------------
#nbelem,nicell,nnodes,nbnodes,npob,npoc,nd
*NBELEM *nelem *npoin *NBNODES *NPOB *nnode *ndime

*# ------------------------------------------
*# izpis koordinat vozlisc
*# ------------------------------------------
# node_number x y
*set elems(all)
*loop nodes
*format {"%8i%15.6e%15.6e%15.6e"}
*NodesNum *NodesCoord(1,real) *NodesCoord(2,real) *NodesCoord(3,real)
*end nodes

*# ------------------------------------------
*# izpis boundary conectivity IBC
*# ------------------------------------------
# Boundary_Element i1 i2 .. i9
*Set Cond boundary *elems *CanRepeat
*loop elems *OnlyInCond
*LoopVar *GlobalNodes(2) *GlobalNodes(5) *GlobalNodes(1) *GlobalNodes(8) *GlobalNodes(4) *GlobalNodes(7) *GlobalNodes(3) *GlobalNodes(6) *GlobalNodes(9)
*end

*# ------------------------------------------
*# izpis domain conectivity IDC
*# ------------------------------------------
# Domain_Cell i1 i2 i3 i4 ... i27
*loop elems
*ElemsNum *ElemsConec
*end elems

*# ------------------------------------------
*# izpis boundary element wall number BEWN
*# ------------------------------------------
# Boundary element wall numbers
# belem,wall
*Set Cond boundary *elems *CanRepeat
*loop elems *OnlyInCond
*LoopVar *cond(1)
*end



