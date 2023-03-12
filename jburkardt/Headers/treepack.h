#ifndef WRAPPER_TREEPACK_H_INCLUDED
#define WRAPPER_TREEPACK_H_INCLUDED

#define vec_next(a,b,c,d) FUNCNAME_VECNEXT(C_DTIPIPB(a,b,c,d))
#define vec_random(a,b,c,d) FUNCNAME_VECRANDOM(C_DTI2PI(a,b,c,d))
#define graph_adj_edge_count(a,b) R_USHRT(FUNCNAME_GRAPHADJEDGECOUNT(C_DTPI(b,a)))
#define graph_adj_is_node_connected(a,b) R_UCHR(FUNCNAME_GRAPHADJISNODECONNECTED(C_DTPI(b,a)))
#define graph_adj_is_tree(a,b) R_UCHR(FUNCNAME_GRAPHADJISTREE(C_DTPI(b,a)))
#define graph_arc_degree(a,b,c,d) FUNCNAME_GRAPHARCDEGREE(C_2DT2PI(a,b,c,d))
#define graph_arc_is_tree(a,b,c,d) R_UCHR(FUNCNAME_GRAPHARCISTREE(C_2DT2PI(a,d,b,c)))
#define graph_arc_node_count(a,b,c) R_USHRT(FUNCNAME_GRAPHARCNODECOUNT(C_DT2PI(a,b,c)))
#define graph_arc_node_max(a,b,c) R_USHRT(FUNCNAME_GRAPHARCNODEMAX(C_DT2PI(a,b,c)))
#define graph_arc_to_graph_adj(a,b,c) FUNCNAME_GRAPHARCTOGRAPHADJ(C_DT2PI(a,b,c))
#define pruefer_to_tree_arc(a,b,c,d) FUNCNAME_PRUEFERTOTREEARC(C_DT3PI(a,b,c,d))
#define pruefer_to_tree_2(a,b,c) FUNCNAME_PRUEFERTOTREE2(C_DT2PI(a,b,c))
#define pruefer_to_tree_2_new(a,b) FUNCNAME_PRUEFERTOTREE2NEW(C_DTPI(a,b))
#define tree_arc_center(a,b,c,d,e,f) FUNCNAME_TREEARCCENTER(C_DT5PI(a,b,c,d,e,f))
#define tree_arc_diam(a,b,c,d,e,f,g) FUNCNAME_TREEARCDIAM(C_DT6PI(a,b,c,d,e,f,g))
#define tree_arc_random(a,b,c,d,e) FUNCNAME_TREEARCRANDOM(C_DT4PI(a,b,c,d,e))
#define tree_arc_to_pruefer(a,b,c) FUNCNAME_TREEARCTOPRUEFER(C_DT2PI(a,b,c))
#define tree_enum(a) R_UINT(FUNCNAME_TREEENUM(C_SUINT(a)))
#define tree_parent_next(a,b,c,d) FUNCNAME_TREEPARENTNEXT(C_DT2PIPB(a,b,c,d))
#define tree_parent_to_arc(a,b,c,d,e) FUNCNAME_TREEPARENTTOARC(C_DT2PIPDTPI(a,b,c,d,e))
#define tree_rb_lex_next(a,b,c) FUNCNAME_TREERBLEXNEXT(C_DTPIPB(a,b,c))
#define tree_rb_to_parent(a,b) FUNCNAME_TREERBTOPARENT(C_DT2PI(a,b))
#define tree_rb_yule(a,b,c) FUNCNAME_TREERBYULE(C_2PIPDT(b,c,a))
#define tree_rooted_code(a,b) FUNCNAME_TREEROOTEDCODE(C_DTPI(a,b))
#define tree_rooted_code_compare(a,b,c,d) R_SHRT(FUNCNAME_TREEROOTEDCODECOMPARE(C_2DT2PI(a,b,c,d)))
#define tree_rooted_depth(a,b,c,d) FUNCNAME_TREEROOTEDDEPTH(C_DT3PI(a,b,c,d))
#define tree_rooted_enum(a) FUNCNAME_TREEROOTEDENUM(C_SUSHRT(a))
#define tree_rooted_random(a,b) FUNCNAME_TREEROOTEDRANDOM(C_DTPI(a,b))

__MATHSUITE __JBURKARDT void * FUNCNAME_VECNEXT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_VECRANDOM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GRAPHADJEDGECOUNT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GRAPHADJISNODECONNECTED(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GRAPHADJISTREE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GRAPHARCISTREE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GRAPHARCNODECOUNT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GRAPHARCNODEMAX(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PRUEFERTOARCTREE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PRUEFERTOTREEARC(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PRUEFERTOTREE2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TREEARCCENTER(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TREEARCDIAM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TREEARCRANDOM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TREEPARENTNEXT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TREEPARENTTOARC(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TREERBLEXNEXT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TREERBYULE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TREEROOTEDCODECOMPARE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TREEROOTEDDEPTH(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GRAPHARCDEGREE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GRAPHARCTOGRAPHADJ(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PRUEFERTOTREE2NEW(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TREEARCTOPRUEFER(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TREERBTOPARENT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TREEROOTEDCODE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TREEROOTEDENUM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TREEROOTEDRANDOM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TREEENUM(void *);

#endif
