#ifndef WRAPPER_GEOMETRY_H_INCLUDED
#define WRAPPER_GEOMETRY_H_INCLUDED

#define dqrank(a,b,c,d,e,f,g,h) FUNCNAME_DQRANK(C_PIT3DTIT2PIPIT(a,b,c,d,e,f,g,h))
#define dqrls(a,b,c,d,e,f,g,h,i,j,k,l) R_SHRT(FUNCNAME_DQRLS(C_PIT3DTITPI3PITPIPITI(a,b,c,d,e,f,g,h,i,j,k,l)))
#define dqrlss(a,b,c,d,e,f,g,h,i,j) FUNCNAME_DQRLSS(C_PIT4DT3PITPIPIT(a,b,c,d,e,f,g,h,i,j))
#define r8mat_amax(a,b,c) R_DBL(FUNCNAME_R8MATAMAX(C_2DTPIT(a,b,c)))
#define qr_solve(a,b,c,d) FUNCNAME_QRSOLVE(C_2DT2PIT(a,b,c,d))
#define r8mat_solve(a,b,c) R_USHRT(FUNCNAME_R8MATSOLVE(C_2DTPIT(a,b,c)))
#define angle_box_2d(a,b,c,d,e,f) FUNCNAME_ANGLEBOX2D(C_IT5PIT(a,b,c,d,e,f))
#define angle_contains_ray_2d(a,b,c,d) R_UCHR(FUNCNAME_ANGLECONTAINSRAY2D(C_PPDBL4(a,b,c,d)))
#define angle_deg_2d(a,b,c) R_DBL(FUNCNAME_ANGLEDEG2D(C_PPDBL3(a,b,c)))
#define angle_half_2d(a,b,c) FUNCNAME_ANGLEHALF2D(C_PPDBL3(a,b,c))
#define angle_rad_2d(a,b,c) R_DBL(FUNCNAME_ANGLERAD2D(C_PPDBL3(a,b,c)))
#define angle_rad_3d(a,b,c) R_DBL(FUNCNAME_ANGLERAD3D(C_PPDBL3(a,b,c)))
#define angle_rad_nd(a,b,c) R_DBL(FUNCNAME_ANGLERADND(C_DT2PIT(a,b,c)))
#define angle_turn_2d(a,b,c) R_DBL(FUNCNAME_ANGLETURN2D(C_PPDBL3(a,b,c)))
#define anglei_deg_2d(a,b,c) R_DBL(FUNCNAME_ANGLEIDEG2D(C_PPDBL3(a,b,c)))
#define anglei_rad_2d(a,b,c) R_DBL(FUNCNAME_ANGLEIRAD2D(C_PPDBL3(a,b,c)))
#define annulus_area_2d(a,b) R_DBL(FUNCNAME_ANNULUSAREA2D(C_PPDBL2(a,b)))
#define annulus_sector_area_2d(a,b,c,d) R_DBL(FUNCNAME_ANNULUSSECTORAREA2D(C_PPDBL4(a,b,c,d)))
#define annulus_sector_centroid_2d(a,b,c,d,e) FUNCNAME_ANNULUSSECTORCENTROID2D(C_4ITPIT(b,c,d,e,a))
#define ball_unit_sample_2d(a) FUNCNAME_BALLUNITSAMPLE2D(a)
#define ball_unit_sample_3d(a) FUNCNAME_BALLUNITSAMPLE3D(a)
#define ball_unit_sample_nd(a,b) FUNCNAME_BALLUNITSAMPLEND(C_DTPI(a,b))
#define basis_map_3d(a,b) R_DBL(FUNCNAME_BASISMAP3D(C_PPDBL2(a,b)))
#define box_01_contains_point_2d(a) R_UCHR(FUNCNAME_BOX01CONTAINSPOINT2D(a))
#define box_01_contains_point_nd(a,b) R_UCHR(FUNCNAME_BOX01CONTAINSPOINTND(C_DTPIT(a,b)))
#define box_contains_point_2d(a,b,c) R_UCHR(FUNCNAME_BOXCONTAINSPOINT2D(C_PPDBL3(a,b,c)))
#define box_contains_point_nd(a,b,c,d) R_UCHR(FUNCNAME_BOXCONTAINSPOINTND(C_DT3PIT(a,b,c,d)))
#define box_ray_int_2d(a,b,c,d,e) FUNCNAME_BOXRAYINT2D(C_PPDBL5(a,b,c,d,e))
#define box_segment_clip_2d(a,b,c,d) R_SHRT(FUNCNAME_BOXSEGMENTCLIP2D(C_PPDBL4(a,b,c,d)))
#define circle_arc_point_near_2d(a,b,c,d,e,f,g) FUNCNAME_CIRCLEARCPOINTNEAR2D(C_3IT4PIT(a,c,d,b,e,f,g))
#define circle_area_2d(a) R_DBL(FUNCNAME_CIRCLEAREA2D(C_SDBL(a)))
#define circle_dia2imp_2d(a,b,c,d) FUNCNAME_CIRCLEDIA2IMP2D(C_PPDBL4(a,b,c,d))
#define circle_exp_contains_point_2d(a,b,c,d) R_USHRT(FUNCNAME_CIRCLEEXPCONTAINSPOINT2D(C_PPDBL4(a,b,c,d)))
#define circle_exp2imp_2d(a,b,c,d,e) FUNCNAME_CIRCLEEXP2IMP2D(C_PPDBL5(a,b,c,d,e))
#define circle_imp_contains_point_2d(a,b,c) R_UCHR(FUNCNAME_CIRCLEIMPCONTAINSPOINT2D(C_DT2PIT(a,b,c)))
#define circle_imp_line_par_int_2d(a,b,c,d,e,f,g,h) FUNCNAME_CIRCLEIMPLINEPARINT2D(C_ITPIT4ITPDTPIT(a,b,c,d,e,f,g,h))
#define circle_imp_point_dist_2d(a,b,c) R_DBL(FUNCNAME_CIRCLEIMPPOINTDIST2D(C_DT2PIT(a,b,c)))
#define circle_imp_point_dist_signed_2d(a,b,c) FUNCNAME_CIRCLEIMPPOINTDISTSIGNED2D(C_DT2PIT(a,b,c))
#define circle_imp_point_near_2d(a,b,c,d) R_DBL(FUNCNAME_CIRCLEIMPPOINTNEAR2D(C_IT3PIT(a,b,c,d)))
#define circle_imp_points_2d(a,b,c) R_DBL(FUNCNAME_CIRCLEIMPPOINTS2D(C_IT3PIT(a,b,c)))
#define circle_imp_points_3d(a,b,c,d) R_DBL(FUNCNAME_CIRCLEIMPPOINTS3D(C_IT2PITDT(a,b,c,d)))
#define circle_imp_points_arc_2d(a,b,c,d,e,f) FUNCNAME_CIRCLEIMPPOINTSARC2D(C_PIT3ITDTPIT(b,a,c,d,e,f))
#define circle_imp2exp_2d(a,b,c,d,e) FUNCNAME_CIRCLEIMP2EXP2D(C_IT4PIT(a,b,c,d,e))
#define circle_llr2imp_2d(a,b,c,d,e) FUNCNAME_CIRCLELLR2IMP2D(C_IT4PIT(e,a,b,c,d))
#define circle_lune_area_2d(a,b,c,d) R_DBL(FUNCNAME_CIRCLELUNEAREA2D(C_ITPIT2IT(a,b,c,d)))
#define circle_lune_centroid_2d(a,b,c,d) FUNCNAME_CIRCLELUNECENTROID2D(C_ITPIT2IT(a,b,c,d))
#define circle_pppr2imp_3d(a,b,c,d,e,f) FUNCNAME_CIRCLEPPPR2IMP3D(C_IT5PIT(d,a,b,c,e,f))
#define circle_ppr2imp_2d(a,b,c) FUNCNAME_CIRCLEPPR2IMP2D(C_IT2PIT(c,a,b))
#define circle_sector_area_2d(a,b,c,d) R_DBL(FUNCNAME_CIRCLESECTORAREA2D(C_ITPIT2IT(a,b,c,d)))
#define circle_sector_centroid_2d(a,b,c,d) FUNCNAME_CIRCLESECTORCENTROID2D(C_ITPIT2IT(a,b,c,d))
#define circle_sector_contains_point_2d(a,b,c,d,e) R_UCHR(FUNCNAME_CIRCLESECTORCONTAINSPOINT2D(C_ITPIT2ITPIT(a,b,c,d,e)))
#define circle_triangle_area_2d(a,b,c,d) R_DBL(FUNCNAME_CIRCLETRIANGLEAREA2D(C_ITPIT2IT(a,b,c,d)))
#define circle_triple_angles_2d(a,b,c,d,e,f) FUNCNAME_CIRCLETRIPLEANGLES2D(C_3IT3PIT(a,b,c,d,e,f))
#define circles_imp_int_2d(a,b,c,d,e,f) FUNCNAME_CIRCLESIMPINT2D(C_2PIT2ITPDTPIT(b,d,c,a,e,f))
#define cone_area_3d(a,b) R_DBL(FUNCNAME_CONEAREA3D(C_PPDBL2(a,b)))
#define cone_centroid_3d(a,b,c) FUNCNAME_CONECENTROID3D(C_IT2PIT(a,b,c))
#define cone_volume_3d(a,b) R_DBL(FUNCNAME_CONEVOLUME3D(C_PPDBL2(a,b)))
#define conv3d(a,b,c,d,e) FUNCNAME_CONV3D(C_CHITDT2PIT(a,b,c,d,e))
#define cot_rad(a) R_DBL(FUNCNAME_COTRAD(C_SDBL(a)))
#define cube_shape_3d(a,b,c,d,e,f) FUNCNAME_CUBESHAPE3D(C_3DTPIT2PI(a,b,c,d,e,f))
#define cube_size_3d(a,b,c,d) FUNCNAME_CUBESIZE3D(C_PPUSHRT4(a,b,c,d))
#define cylinder_point_dist_3d(a,b,c,d) R_DBL(FUNCNAME_CYLINDERPOINTDIST3D(C_IT3PIT(c,a,b,d)))
#define cylinder_point_dist_signed_3d(a,b,c,d) R_DBL(FUNCNAME_CYLINDERPOINTDISTSIGNED3D(C_IT3PIT(c,a,b,d)))
#define cylinder_point_inside_3d(a,b,c,d) R_UCHR(FUNCNAME_CYLINDERPOINTINSIDE3D(C_IT3PIT(c,a,b,d)))
#define cylinder_point_near_3d(a,b,c,d) FUNCNAME_CYLINDERPOINTNEAR3D(C_IT3PIT(c,a,b,d))
#define cylinder_sample_3d(a,b,c,d,e) FUNCNAME_CYLINDERSAMPLE3D(C_2PITITDTPI(a,b,c,d,e))
#define cylinder_volume_3d(a,b,c) R_DBL(FUNCNAME_CYLINDERVOLUME3D(CC_IT2PIT(c,a,b)))
#define degrees_to_radians(a) R_DBL(FUNCNAME_DEGREESTORADIANS(C_SDBL(a)))
#define dge_det(a,b,c) R_DBL(FUNCNAME_DGEDET(C_DTPITPDT(a,b,c)))
#define dge_fa(a,b,c) R_USHRT(FUNCNAME_DGE0FA(C_DTPITPDT(a,b,c)))
#define dge_sl(a,b,c,d,e) FUNCNAME_DGE0SL(C_PIT2DTPDTPIT(b,d,e,a,c))
#define direction_pert_3d(a,b,c) FUNCNAME_DIRECTIONPERT3D(C_ITPITPI(a,b,c))
#define direction_uniform_2d(a) FUNCNAME_DIRECTIONUNIFORM2D(a)
#define direction_uniform_3d(a) FUNCNAME_DIRECTIONUNIFORM3D(a)
#define direction_uniform_nd(a,b) FUNCNAME_DIRECTIONUNIFORMND(C_DTPI(a,b))
#define disk_point_dist_3d(a,b,c,d) R_DBL(FUNCNAME_DISKPOINTDIST3D(C_IT3PIT(b,a,c,d)))
#define dms_to_radians(a,b,c) R_DBL(FUNCNAME_DMSTORADIANS(C_PUSHRT3(a,b,c)))
#define dodec_shape_3d(a,b,c,d,e,f) FUNCNAME_DODECSHAPE3D(C_3DTPIT2PI(a,b,c,d,e,f))
#define dual_shape_3d(a,b,c,d,e,f,g,h,i,j,k,l) FUNCNAME_DUALSHAPE3D(C_3DTPITPDTPI3DTPIT2PDT(a,b,c,d,e,f,g,h,i,j,k,l))
#define dual_size_3d(a,b,c,d,e,f,g,h,i,j) FUNCNAME_DUALSIZE3D(C_3DTPIT6PDT(a,b,c,d,e,f,g,h,i,j))
#define ellipse_area_2d(a,b) R_DBL(FUNCNAME_ELLIPSEAREA2D(C_PDBL2(a,b)))
#define ellipse_point_dist_2d(a,b,c) R_DBL(FUNCNAME_ELLIPSEPOINTDIST2D(C_2ITPIT(a,b,c)))
#define ellipse_point_near_2d(a,b,c) FUNCNAME_ELLIPSEPOINTNEAR2D(C_2ITPIT(a,b,c))
#define ellipse_points_2d(a,b,c,d,e,f) FUNCNAME_ELLIPSEPOINTS2D(C_PIT3ITDTPIT(a,b,c,d,e,f))
#define ellipse_points_arc_2d(a,b,c,d,e,f) FUNCNAME_ELLIPSEPOINTSARC2D(C_PIT5ITDTPIT(a,b,c,d,e,f))
#define enorm0_nd(a,b,c) R_DBL(FUNCNAME_ENORM0ND(C_ITPITPI(a,b,c)))
#define get_seed() R_INT(FUNCNAME_GETSEED(NULL))
#define glob2loc_3d(a,b,c,d,e,f,g,h,i) FUNCNAME_GLOB2LOC3D(C_6IT3PIT(a,b,c,d,e,f,g,h,i))
#define halfplane_contains_point_2d(a,b,c) R_UCHR(FUNCNAME_HALFPLANECONTAINSPOINT2D(C_PPDBL3(a,b,c)))
#define halfspace_imp_triangle_int_3d(a,b,c,d,e,f) R_USHRT(FUNCNAME_HALFSPACEIMPTTRIANGLEINT3D(C_4IT2PIT(a,b,c,d,e,f)))
#define halfspace_norm_triangle_int_3d(a,b,c,d) R_USHRT(FUNCNAME_HALFSPACENORMTRIANGLEINT3D(C_PPDBL4(a,b,c,d)))
#define halfspace_triangle_int_3d(a,b,c,d,e) R_USHRT(FUNCNAME_HALFSPACETRIANGLEINT3D(C_ITPIT2ITPIT(a,d,b,c,e)))
#define haversine(a) R_DBL(FUNCNAME_HAVERSINE(C_SDBL(a)))
#define helix_shape_3d(a,b,c,d,e,f) FUNCNAME_HELIXSHAPE3D(C_ITDT3ITPIT(a,b,c,d,e,f))
#define hexagon_area_2d(a) R_DBL(FUNCNAME_HEXAGONAREA2D(C_SDBL(a)))
#define hexagon_contains_point_2d(a,b) R_UCHR(FUNCNAME_HEXAGONCONTAINSPOINT2D(C_PDBL2(a,b)))
#define hexagon_shape_2d(a,b) FUNCNAME_HEXAGONSHAPE2D(C_ITPIT(a,b))
#define hexagon_unit_area_2d() R_DBL(FUNCNAME_HEXAGONUNITAREA2D(NULL))
#define hexagon_vertices_2d(a) FUNCNAME_HEXAGONVERTICES2D(a)
#define i4_dedekind_factor(a,b) R_DBL(FUNCNAME_I4DEDEKINDFACTOR(C_PUSHRT2(a,b)))
#define i4_dedekind_sum(a,b) R_DBL(FUNCNAME_I4DEDEKINDSUM(C_PUSHRT2(a,b)))
#define i4_factorial(a) R_INT(FUNCNAME_I4FACTORIAL(C_SINT(a)))
#define i4_factorial2(a) R_INT(FUNCNAME_I4FACTORIAL2(C_SINT(a)))
#define i4_gcd(a,b) R_INT(FUNCNAME_I4GCD(C_PINT2(a,b)))
#define i4_modp(a,b) R_INT(FUNCNAME_I4MODP(C_PINT2(a,b)))
#define i4_sign(a) R_INT(FUNCNAME_I4SIGN(C_SINT(a)))
#define i4_swap(a,b) FUNCNAME_I4SWAP(C_PPINT2(a,b))
#define icos_shape(a,b,c,d,e,f,g,h) FUNCNAME_ICOSSHAPE(C_4DTPIT3PI(a,b,c,d,e,f,g,h))
#define icos_size(a,b,c,d) FUNCNAME_ICOSSIZE(C_PPUSHRT4(a,b,c,d))
#define line_exp_is_degenerate_nd(a,b,c) R_UCHR(FUNCNAME_LINEEXPISDEGENERATEND(C_DT2PIT(a,b,c)))
#define line_exp_normal_2d(a,b) FUNCNAME_LINEEXPNORMAL2D(C_PPDBL2(a,b))
#define line_exp_perp_2d(a,b,c,d) FUNCNAME_LINEEXPPERP2D(C_3PITPB(a,b,c,d))
#define line_exp_point_dist_2d(a,b,c) R_DBL(FUNCNAME_LINEEXPPOINTDIST2D(C_PPDBL3(a,b,c)))
#define line_exp_point_dist_3d(a,b,c) R_DBL(FUNCNAME_LINEEXPPOINTDIST3D(C_PPDBL3(a,b,c)))
#define line_exp_point_dist_signed_2d(a,b,c) R_DBL(FUNCNAME_LINEEXPPOINTDISTSIGNED2D(C_PPDBL3(a,b,c)))
#define line_exp_point_near_2d(a,b,c,d,e,f) FUNCNAME_LINEEXPPOINTNEAR2D(C_PPDBL6(a,b,c,d,e,f))
#define line_exp_point_near_3d(a,b,c,d,e,f) FUNCNAME_LINEEXPPOINTNEAR3D(C_PPDBL6(a,b,c,d,e,f))
#define line_exp2imp_2d(a,b,c,d,e) FUNCNAME_LINEEXP2IMP2D(C_PPDBL5(a,b,c,d,e))
#define line_exp2par_2d(a,b,c,d,e,f) FUNCNAME_LINEEXP2PAR2D(C_PPDBL6(a,b,c,d,e,f))
#define line_exp2par_3d(a,b,c,d,e,f,g,h) FUNCNAME_LINEEXP2PAR3D(C_PPDBL8(a,b,c,d,e,f,g,h))
#define line_imp_is_degenerate_2d(a,b,c) R_UCHR(FUNCNAME_LINEIMPISDEGENERATE2D(C_PDBL3(a,b,c)))
#define line_imp_point_dist_2d(a,b,c,d) R_DBL(FUNCNAME_LINEIMPPOINTDIST2D(C_ITPIT2IT(a,d,b,c)))
#define line_imp_point_dist_signed_2d(a,b,c,d) R_DBL(FUNCNAME_LINEIMPPOINTDISTSIGNED2D(C_ITPIT2IT(b,a,c,d)))
#define line_imp2exp_2d(a,b,c,d,e) FUNCNAME_LINEIMP2EXP2D(_ITPIT2ITPIT(a,d,b,c,e))
#define line_imp2par_2d(a,b,c,d,e,f,g) FUNCNAME_LINEIMP2PAR2D(C_3IT4PIT(a,b,c,d,e,f,g))
#define line_par_point_dist_2d(a,b,c,d,e) R_DBL(FUNCNAME_LINEPARPOINTDIST2D(C_4ITPIT(a,b,c,d,e)))
#define line_par_point_dist_3d(a,b,c,d,e,f,g) R_DBL(FUNCNAME_LINEPARPOINTDIST3D(C_6ITPIT(a,b,c,d,e,f,g)))
#define line_par_point_near_2d(a,b,c,d,e) FUNCNAME_LINEPARPOINTNEAR2D(C_4ITPIT(a,b,c,d,e))
#define line_par_point_near_3d(a,b,c,d,e,f,g) FUNCNAME_LINEPARPOINTNEAR3D(C_6ITPIT(a,b,c,d,e,f,g))
#define line_par2exp_2d(a,b,c,d,e,f) FUNCNAME_LINEPAR2EXP2D(C_4IT2PIT(a,b,c,d,e,f))
#define line_par2exp_3d(a,b,c,d,e,f,g,h) FUNCNAME_LINEPAR2EXP3D(C_6IT2PIT(a,b,c,d,e,f,g,h))
#define line_par2imp_2d(a,b,c,d,e,f,g) FUNCNAME_LINEPAR2IMP2D(C_4IT3PIT(a,b,c,d,e,f,g))
#define lines_exp_angle_3d(a,b,c,d) R_DBL(FUNCNAME_LINESEXPANGLE3D(C_PPDBL4(a,b,c,d)))
#define lines_exp_angle_nd(a,b,c,d,e) R_DBL(FUNCNAME_LINESEXPANGLEND(C_4ITPIT(a,b,c,d,e)))
#define lines_exp_dist_3d(a,b,c,d) R_DBL(FUNCNAME_LINESEXPDIST3D(C_PPDBL4(a,b,c,d)))
#define lines_exp_dist_3d_2(a,b,c,d) R_DBL(FUNCNAME_LINESEXPDIST3D2(C_PPDBL4(a,b,c,d)))
#define lines_exp_equal_2d(a,b,c,d) R_UCHR(FUNCNAME_LINESEXPEQUAL2D(C_PPDBL4(a,b,c,d)))
#define lines_exp_int_2d(a,b,c,d,e,f) FUNCNAME_LINESEXPINT2D(C_4PITPDTPIT(a,b,c,d,e,f))
#define lines_exp_near_3d(a,b,c,d,e,f) FUNCNAME_LINESEXPNEAR3D(C_PPDBL6(a,b,c,d,e,f))
#define lines_exp_parallel_2d(a,b,c,d) R_UCHR(FUNCNAME_LINESEXPPARALLEL2D(C_PPDBL4(a,b,c,d)))
#define lines_exp_parallel_3d(a,b,c,d) R_UCHR(FUNCNAME_LINESEXPPARALLEL3D(C_PPDBL4(a,b,c,d)))
#define lines_imp_angle_2d(a,b,c,d,e,f) R_DBL(FUNCNAME_LINESIMPANGLE2D(C_PDBL6(a,b,c,d,e,f)))
#define lines_imp_dist_2d(a,b,c,d,e,f) R_DBL(FUNCNAME_LINESIMPDIST2D(C_PDBL6(a,b,c,d,e,f)))
#define lines_imp_int_2d(a,b,c,d,e,f,g,h) FUNCNAME_LINESIMPINT2D(C_6ITPDTPIT(a,b,c,d,e,f,g,h))
#define lines_par_angle_2d(a,b,c,d,e,f,g,h) R_DBL(FUNCNAME_LINESPARANGLE3D(C_PDBL8(a,b,c,d,e,f,g,h)))
#define lines_par_dist_3d(a,b,c,d,e,f,g,h,i,j,k,l) R_DBL(FUNCNAME_LINESPARDIST3D(C_PDBL12(a,b,c,d,e,f,g,h,i,j,k,l)))
#define lines_par_int_2d(a,b,c,d,e,f,g,h,i,j,k) FUNCNAME_LINESPARINT2D(C_8IT3PIT(a,b,c,d,e,f,g,h,i,j,k))
#define loc2glob_3d(a,b,c,d,e,f,g,h,i) FUNCNAME_LOC2GLOB3D(C_6IT3PIT(a,b,c,d,e,f,g,h,i))
#define minabs(a,b,c,d,e,f,g,h) FUNCNAME_MINABS(C_6IT2PIT(a,b,c,d,e,f,g,h))
#define minquad(a,b,c,d,e,f,g,h) R_USHRT(FUNCNAME_MINQUAD(C_6IT2PIT(a,b,c,d,e,f,g,h)))
#define octahedron_shape_3d(a,b,c,d,e,f) FUNCNAME_OCTAHEDRONSHAPE3D(C_3DTPIT2PI(a,b,c,d,e,f))
#define octahedron_size_3d(a,b,c,d) FUNCNAME_OCTAHEDRONSIZE3D(C_PPUSHRT4(a,b,c,d))
#define parabola_ex(a,b,c,d,e,f,g,h) R_USHRT(FUNCNAME_PARABOLAEX(C_6IT2PIT(a,b,c,d,e,f,g,h)))
#define parabola_ex2(a,b,c,d,e,f,g,h,i,j,k) R_USHRT(FUNCNAME_PARABOLAEX2(C_6IT5PIT(a,b,c,d,e,f,g,h,i,j,k)))
#define parallelogram_area_2d(a) R_DBL(FUNCNAME_PARALLELOGRAMAREA2D(a))
#define parallelogram_area_3d(a) R_DBL(FUNCNAME_PARALLELOGRAMAREA3D(a))
#define parallelogram_contains_point_2d(a,b,c,d) R_UCHR(FUNCNAME_PARALLELOGRAMCONTAINSPOINT2D(C_PPDBL4(a,b,c,d)))
#define parallelogram_contains_point_3d(a,b,c,d) R_UCHR(FUNCNAME_PARALLELOGRAMCONTAINSPOINT3D(C_PPDBL4(a,b,c,d)))
#define parallelogram_point_dist_3d(a,b,c,d) R_DBL(FUNCNAME_PARALLELOGRAMPOINTDIST3D(C_PPDBL4(a,b,c,d)))
#define parallelepiped_contains_point_3d(a,b,c,d,e) R_UCHR(FUNCNAME_PARALLEPIPEDCONTAINSPOINT3D(C_PPDBL5(a,b,c,d,e)))
#define parallelepiped_point_dist_3d(a,b,c,d,e) R_DBL(FUNCNAME_PARALLEPIPEDPOINTDIST3D(C_PPDBL5(a,b,c,d,e)))
#define plane_exp_point_dist_3d(a,b,c,d,e) R_DBL(FUNCNAME_PLANEEXPPOINTDIST3D(C_PPDBL5(a,b,c,d,e)))
#define plane_exp_normal_3d(a,b,c,d) FUNCNAME_PLANEEXPNORMAL3D(C_PPDBL4(a,b,c,d))
#define plane_exp_pro2(a,b,c,d,e,f,g) FUNCNAME_PLANEEXPPRO2(C_3PITDT3PIT(a,b,c,d,e,f,g))
#define plane_exp_pro3(a,b,c,d,e,f) FUNCNAME_PLANEEXPPRO3(C_DT5PIT(d,a,b,c,e,f))
#define plane_exp_project_3d(a,b,c,d,e,f,g,h) FUNCNAME_PLANEEXPPROJECT3D(C_4PITDT2PITPDT(a,b,c,d,e,f,g,h))
#define plane_exp2imp_3d(a,b,c,d,e,f,g) FUNCNAME_PLANEEXP2IMP3D(C_PPDBL7(a,b,c,d,e,f,g))
#define plane_exp2normal_3d(a,b,c,d,e) FUNCNAME_PLANEEXP2NORMAL3D(C_PPDBL5(a,b,c,d,e))
#define plane_imp_is_degenerate_3d(a,b,c) R_UCHR(FUNCNAME_PLANEIMPISDEGENERATE3D(C_PDBL3(a,b,c)))
#define plane_imp_line_par_int_3d(a,b,c,d,e,f,g,h,i,j,k) R_UCHR(FUNCNAME_PLANEIMPLINEPARINT3D(C_10ITPIT(a,b,c,d,e,f,g,h,i,j,k)))
#define plane_imp_point_dist_3d(a,b,c,d,e) R_DBL(FUNCNAME_PLANEIMPPOINTDIST3D(C_4ITPIT(a,b,c,d,e)))
#define plane_imp_point_dist_signed_3d(a,b,c,d,e) R_DBL(FUNCNAME_PLANEIMPPOINTDISTSIGNED3D(C_4ITPIT(a,b,c,d,e)))
#define plane_imp_point_near_3d(a,b,c,d,e,f) FUNCNAME_PLANEIMPPOINTNEAR3D(C_4IT2PIT(a,b,c,d,e,f))
#define plane_imp_segment_near_3d(a,b,c,d,e,f,g,h,i) FUNCNAME_PLANEIMPSEGMENTNEAR3D(C_2PIT4IT3PIT(a,b,c,d,e,f,g,h,i))
#define plane_imp_triangle_int_3d(a,b,c,d,e,f,g) FUNCNAME_PLANEIMPTRIANGLEINT3D(C_4ITPITPIPIT(a,b,c,d,e,f,g))
#define plane_imp_triangle_int_add_3d(a,b,c,d,e,f) FUNCNAME_PLANEIMPTRIANGLEINTADD3D(C_2PIT2ITPDTPIT(a,b,c,d,e,f))
#define plane_imp_triangle_near_3d(a,b,c,d,e,f,g) R_USHRT(FUNCNAME_PLANEIMPTRIANGLENEAR3D(C_4IT3PIT(b,c,d,e,f,g,a)))
#define plane_imp2exp_3d(a,b,c,d,e,f,g) FUNCNAME_PLANEIMP2EXP3D(C_4IT3PIT(a,b,c,d,e,f,g))
#define plane_imp2normal_3d(a,b,c,d,e,f) FUNCNAME_PLANEIMP2NORMAL3D(C_4IT2PIT(a,b,c,d,e,f))
#define plane_normal_basis_3d(a,b,c,d) FUNCNAME_PLANENORMALBASIS3D(C_PPDBL4(a,b,c,d))
#define plane_normal_line_exp_int_3d(a,b,c,d,e) R_INT(FUNCNAME_PLANENORMALLINEEXPINT3D(C_PPDBL5(a,b,c,d,e)))
#define plane_normal_qr_to_xyz(a,b,c,d,e,f) FUNCNAME_PLANENORMALQRTOXYZ(C_DT5PIT(e,a,b,c,d,f))
#define plane_normal_tetrahedron_intersect(a,b,c,d,e) FUNCNAME_PLANENORMALTETRAHEDRONINTERSECT(C_PDT4PIT(d,a,b,c,e))
#define plane_normal_triangle_int_3d(a,b,c,d) R_USHRT(FUNCNAME_PLANENORMALTRIANGLEINT3D(C_PPDBL4(a,b,c,d)))
#define plane_normal_uniform_3d(a,b,c) FUNCNAME_PLANENORMALUNIFORM3D(C_PI2PIT(a,b,c))
#define plane_normal_uniform_nd(d,a,b,c) FUNCNAME_PLANENORMALUNIFORMND(C_PITDTPIPIT(d,a,b,c))
#define plane_normal2exp_3d(a,b,c,d,e) FUNCNAME_PLANENORMAL2EXP3D(C_PPDBL5(a,b,c,d,e))
#define plane_normal2imp_3d(a,b,c,d,e,f) FUNCNAME_PLANENORMAL2IMP3D(C_PPDBL6(a,b,c,d,e,f))
#define planes_imp_angle_3d(a,b,c,d,e,f,g,h) R_DBL(FUNCNAME_PLANESIMPANGLE3D(C_PDBL8(a,b,c,d,e,f,g,h)))
#define points_avoid_point_naive_2d(a,b,c) R_UCHR(FUNCNAME_POINTSAVOIDPOINTNAIVE2D(C_DT2PIT(a,b,c)))
#define points_bisect_line_imp_2d(a,b,c,d,e) FUNCNAME_POINTSBISECTLINEIMP2D(C_PPDBL5(a,b,c,d,e))
#define points_bisect_line_par_2d(a,b,c,d,e,f) FUNCNAME_POINTSBISECTLINEPAR2D(C_PDBL6(a,b,c,d,e,f))
#define points_centroid_2d(a,b) R_USHRT(FUNCNAME_POINTSCENTROID2D(C_DTPIT(a,b)))
#define points_colin_2d(a,b,c) R_DBL(FUNCNAME_POINTSCOLIN2D(C_PPDBL3(a,b,c)))
#define points_colin_3d(a,b,c) R_DBL(FUNCNAME_POINTSCOLIN3D(C_PPDBL3(a,b,c)))
#define points_dist_2d(a,b) R_DBL(FUNCNAME_POINTSDIST2D(C_PPDBL2(a,b)))
#define points_dist_3d(a,b) R_DBL(FUNCNAME_POINTSDIST3D(C_PPDBL2(a,b)))
#define points_dist_nd(a,b,c) R_DBL(FUNCNAME_POINTSDISTND(C_DT2PIT(a,b,c)))
#define points_hull_2d(a,b,c,d) FUNCNAME_POINTSHULL2D(C_DTPITPDTPI(a,b,c,d))
#define points_point_near_naive_2d(a,b,c,d) R_USHRT(FUNCNAME_POINTSPOINTNEARNAIVE2D(C_DT3PIT(a,b,c,d)))
#define points_point_near_naive_3d(a,b,c,d) R_USHRT(FUNCNAME_POINTSPOINTNEARNAIVE3D(C_DT3PIT(a,b,c,d)))
#define points_point_near_naive_nd(a,b,c,d,e) R_USHRT(FUNCNAME_POINTSPOINTNEARNAIVEND(C_2DT3PIT(a,b,c,d,e)))
#define points_points_near_naive_2d(a,b,c,d) FUNCNAME_POINTSPOINTSNEARNAIVE2D(C_2DT2PIT(a,c,b,d))
#define points_points_near_naive_3d(a,b,c,d) FUNCNAME_POINTSPOINTSNEARNAIVE3D(C_2DT2PIT(a,c,b,d))
#define polar_to_xy(a,b,c) FUNCNAME_POLARTOXY(C_2ITPIT(a,b,c))
#define polygon_1_2d(a,b) R_DBL(FUNCNAME_POLYGON12D(C_DTPIT(a,b)))
#define polygon_angles_2d(a,b) FUNCNAME_POLYGONANGLES2D(C_DTPIT(a,b))
#define polygon_area_2d(a,b) R_DBL(FUNCNAME_POLYGONAREA2D(C_DTPIT(a,b)))
#define polygon_area_2d_2(a,b) R_DBL(FUNCNAME_POLYGONAREA2D2(C_DTPIT(a,b)))
#define polygon_area_3d(a,b,c) R_DBL(FUNCNAME_POLYGONAREA3D(C_DT2PIT(a,b,c)))
#define polygon_area_3d_2(a,b,c) R_DBL(FUNCNAME_POLYGONAREA3D2(C_DT2PIT(a,b,c)))
#define polygon_centroid_2d(a,b) FUNCNAME_POLYGONCENTROID2D(C_DTPIT(a,b))
#define polygon_centroid_2d_2(a,b) FUNCNAME_POLYGONCENTROID2D2(C_DTPIT(a,b))
#define polygon_centroid_3d(a,b) FUNCNAME_POLYGONCENTROID3D(C_DTPIT(a,b))
#define polygon_contains_point_2d(a,b,c) R_UCHR(FUNCNAME_POLYGONCONTAINSPOINT2D(C_DT2PIT(a,b,c)))
#define polygon_contains_point_2d_2(a,b,c) R_UCHR(FUNCNAME_POLYGONCONTAINSPOINT2D2(C_DT2PIT(a,b,c)))
#define polygon_diameter_2d(a,b) R_DBL(FUNCNAME_POLYGONDIAMETER2D(C_DTPIT(a,b)))
#define polygon_expand_2d(a,b,c) FUNCNAME_POLYGONEXPAND2D(C_DTPITIT(a,b,c))
#define polygon_inrad_data_2d(a,b,c,d,e) FUNCNAME_POLYGONINRADDATA2D(C_DTIT3PIT(a,b,c,d,e))
#define polygon_is_convex(a,b) R_USHRT(FUNCNAME_POLYGONISCONVEX(C_DTPIT(a,b)))
#define polygon_lattice_area_2d(a,b) R_DBL(FUNCNAME_POLYGONLATTICEAREA2D(C_PUSHRT2(a,b)))
#define polygon_normal_3d(a,b) FUNCNAME_POLYGONNORMAL3D(C_DTPIT(a,b))
#define polygon_outrad_data_2d(a,b,c,d,e) FUNCNAME_POLYGONOUTRADDATA2D(C_DTIT3PIT(a,b,c,d,e))
#define polygon_side_data_2d(a,b,c,d,e) FUNCNAME_POLYGONSIDEDATA2D(C_DTIT3PIT(a,b,c,d,e))
#define polygon_solid_angle_3d(a,b,c) R_DBL(FUNCNAME_POLYGONSOLIDANGLE3D(C_DT2PIT(a,b,c)))
#define polygon_x_2d(a,b) R_DBL(FUNCNAME_POLYGONX2D(C_DTPIT(a,b)))
#define polygon_y_2d(a,b) R_DBL(FUNCNAME_POLYGONY2D(C_DTPIT(a,b)))
#define polygon_xx_2d(a,b) R_DBL(FUNCNAME_POLYGONXX2D(C_DTPIT(a,b)))
#define polygon_xy_2d(a,b) R_DBL(FUNCNAME_POLYGONXY2D(C_DTPIT(a,b)))
#define polygon_yy_2d(a,b) R_DBL(FUNCNAME_POLYGONYY2D(C_DTPIT(a,b)))
#define polyhedron_area_3d(a,b,c,d,e,f) R_DBL(FUNCNAME_POLYHEDRONAREA3D(C_PIT2DTPDTDTPDT(a,b,c,d,e,f)))
#define polyhedron_centroid_3d(a,b,c,d,e,f) FUNCNAME_POLYHEDRONCENTROID3D(C_PIT2DTPDTDTPDT(a,b,c,d,e,f))
#define polyhedron_contains_point_3d(a,b,c,d,e,f,g) R_UCHR(FUNCNAME_POLYHEDRONCONTAINSPOINT3D(C_PIT2DTPDTDTPDT(a,b,c,d,e,f,g)))
#define polyhedron_volume_3d(a,b,c,d,e,f) R_DBL(FUNCNAME_POLYHEDRONVOLUME3D(C_PIT2DTPDTDTPDT(a,b,c,d,e,f)))
#define polyhedron_volume_3d_2(a,b,c,d,e,f) R_DBL(FUNCNAME_POLYHEDRONVOLUME3D2(C_PIT2DTPDTDTPDT(a,b,c,d,e,f)))
#define polyline_arclength_nd(a,b,c) FUNCNAME_POLYLINEARCLENGTHND(C_2DTPIT(a,b,c))
#define polyline_index_point_nd(a,b,c,d) FUNCNAME_POLYLINEINDEXPOINTND(C_2DTPITIT(a,b,c,d))
#define polyline_length_nd(a,b,c) R_DBL(FUNCNAME_POLYLINELENGTHND(C_2DTPIT(a,b,c)))
#define polyline_points_nd(a,b,c,d) FUNCNAME_POLYLINEPOINTSND(C_2DTPITIT(a,b,c,d))
#define polyloop_arclength_nd(a,b,c) FUNCNAME_POLYLOOPARCLENGTHND(C_2DTPIT(a,b,c))
#define polyloop_points_nd(a,b,c,d) FUNCNAME_POLYLOOPPOINTSND(C_2DTPITIT(a,b,c,d))
#define provec(a,b,c,d,e,f) FUNCNAME_PROVEC(C_2DT4PIT(a,b,c,d,e,f))
#define pyramid_volume_3d(a,b) R_DBL(FUNCNAME_PYRAMIDVOLUME3D(C_PDBL2(a,b)))
#define quad_area_2d(a) R_DBL(FUNCNAME_QUADAREA2D(a))
#define quad_area2_2d(a) R_DBL(FUNCNAME_QUADAREA2D2(a))
#define quad_area_3d(a) R_DBL(FUNCNAME_QUADAREA3D(a))
#define quad_contains_point_2d(a,b) R_UCHR(FUNCNAME_QUADCONTAINSPOINT2D(C_PPDBL2(a,b)))
#define quad_point_dist_2d(a,b) R_DBL(FUNCNAME_QUADPOINTDIST2D(C_PPDBL2(a,b)))
#define quad_point_dist_signed_2d(a,b) R_DBL(FUNCNAME_QUADPOINTDISTSIGNED2D(C_PPDBL2(a,b)))
#define quad_point_near_2d(a,b,c,d) FUNCNAME_QUADPOINTNEAR2D(C_PPDBL4(a,b,c,d))
#define quat_conj(a) FUNCNAME_QUATCONJ(a)
#define quat_inv(a) FUNCNAME_QUATINV(a)
#define quat_mul(a,b) FUNCNAME_QUATMUL(C_PPDBL2(a,b))
#define quat_norm(a) R_DBL(FUNCNAME_QUATNORM(a))
#define r8mat_mv_new(a,b,c,d) FUNCNAME_R8MATMVNEW(&((_2dt2pit){(a),(b),(c),(d)}))
#define r8mat_solve_2d(a,b,c) FUNCNAME_R8MATSOLVE2D(C_PPDBL3(a,b,c))
#define r8vec_angle_3d(a,b) R_DBL(FUNCNAME_R8VECANGLE3D(C_PPDBL2(a,b)))
#define r8vec_bracket(a,b,c,d,e) FUNCNAME_R8VECBRACKET(C_DTPITIT2PDT(a,b,c,d,e))
#define radec_distance_3d(a,b,c,d) R_DBL(FUNCNAME_RADECDISTANCE3D(C_PDBL4(a,b,c,d)))
#define radec_to_xyz(a,b) FUNCNAME_RADECTOXYZ(C_PDBL2(a,b))
#define radians_to_degrees(a) R_DBL(FUNCNAME_RADIANSTODEGREES(C_SDBL(a)))
#define radians_to_dms(a,b,c,d) FUNCNAME_RADIANSTODMS(C_IT3PI(a,b,c,d))
#define rotation_axis_vector_3d(a,b,c,d) FUNCNAME_ROTATIONAXISVECTOR3D(C_IT3PIT(b,a,c,d))
#define rotation_axis2mat_3d(a,b,c) FUNCNAME_ROTATIONAXIS2MAT3D(C_IT2PIT(b,a,c))
#define rotation_mat_vector_3d(a,b,c) FUNCNAME_ROTATIONMATVECTOR3D(C_PPDBL3(a,b,c))
#define rotation_mat2axis_3d(a,b,c) FUNCNAME_ROTATIONMAT2AXIS3D(C_PPDBL3(a,b,c))
#define rotation_mat2quat_3d(a,b) FUNCNAME_ROTATIONMAT2QUAT3D(C_PPDBL2(a,b))
#define rotation_quat_vector_3d(a,b,c) FUNCNAME_ROTATIONQUATVECTOR3D(C_PPDBL3(a,b,c))
#define rotation_quat2axis_3d(a,b,c) FUNCNAME_ROTATIONQUAT2AXIS3D(C_PPDBL3(a,b,c))
#define rotation_quat2mat_3d(a,b) FUNCNAME_ROTATIONQUAT2MAT3D(C_PPDBL2(a,b))
#define rtp_to_xyz(a,b,c,d) FUNCNAME_RTPTOXYZ(C_ITPIT2IT(a,d,b,c))
#define segment_contains_point_1d(a,b,c,d) FUNCNAME_SEGMENTCONTAINSPOINT1D(C_ITPIT2IT(a,d,b,c))
#define segment_contains_point_2d(a,b,c,d) FUNCNAME_SEGMENTCONTAINSPOINT2D(C_PPDBL4(a,b,c,d))
#define segment_point_coords_2d(a,b,c,d,e) FUNCNAME_SEGMENTPOINTCOORDS2D(C_PPDBL5(a,b,c,d,e))
#define segment_point_coords_3d(a,b,c,d,e) FUNCNAME_SEGMENTPOINTCOORDS3D(C_PPDBL5(a,b,c,d,e))
#define segment_point_dist_2d(a,b,c) R_DBL(FUNCNAME_SEGMENTPOINTDIST2D(C_PPDBL3(a,b,c)))
#define segment_point_dist_3d(a,b,c) R_DBL(FUNCNAME_SEGMENTPOINTDIST3D(C_PPDBL3(a,b,c)))
#define segment_point_near_2d(a,b,c,d,e,f) FUNCNAME_SEGMENTPOINTNEAR2D(C_PPDBL6(a,b,c,d,e,f))
#define segment_point_near_3d(a,b,c,d,e,f) FUNCNAME_SEGMENTPOINTNEAR3D(C_PPDBL6(a,b,c,d,e,f))
#define segments_curvature_2d(a,b,c) R_DBL(FUNCNAME_SEGMENTSCURVATURE2D(C_PPDBL3(a,b,c)))
#define segments_dist_2d(a,b,c,d) R_DBL(FUNCNAME_SEGMENTSDIST2D(C_PPDBL4(a,b,c,d)))
#define segments_dist_3d(a,b,c,d) R_DBL(FUNCNAME_SEGMENTSDIST3D(C_PPDBL4(a,b,c,d)))
#define segments_dist_3d_old(a,b,c,d) R_DBL(FUNCNAME_SEGMENTSDIST3DOLD(C_PPDBL4(a,b,c,d)))
#define segments_int_1d(a,b,c,d,e,f) R_DBL(FUNCNAME_SEGMENTSINT1D(C_4IT2PIT(a,b,c,d,e,f)))
#define segments_int_2d(a,b,c,d,e,f) R_DBL(FUNCNAME_SEGMENTSINT2D(C_4PITPDTPIT(a,b,c,d,e,f)))
#define shape_point_dist_2d(a,b,c,d) R_DBL(FUNCNAME_SHAPEPOINTDIST2D(C_DT3PIT(c,d,a,b)))
#define shape_point_near_2d(a,b,c,d,e,f) FUNCNAME_SHAPEPOINTNEAR2D(C_DT5PIT(c,a,b,d,e,f))
#define shape_ray_int_2d(a,b,c,d,e,f) FUNCNAME_SHAPERAYINT2D(C_DT5PIT(c,a,b,d,e,f))
#define simplex_lattice_layer_point_next(a,b,c,d) FUNCNAME_SIMPLEXLATTICELAYERPOINTNEXT(C_DTPIPDTPB(a,b,c,d))
#define simplex_lattice_point_next(a,b,c,d) FUNCNAME_SIMPLEXLATTICEPOINTNEXT(C_DT2PIPB(a,b,c,d))
#define simplex_unit_lattice_point_nd(a,b) R_USHRT(FUNCNAME_SIMPLEXUNITLATTICEPOINTND(C_PUSHRT2(a,b)))
#define simplex_unit_volume_nd(a) R_DBL(FUNCNAME_SIMPLEXUNITVOLUMEND(C_SUSHRT(a)))
#define simplex_volume_nd(a,b) R_DBL(FUNCNAME_SIMPLEXVOLUMEND(C_DTPIT(a,b)))
#define soccer_size_3d(a,b,c,d) FUNCNAME_SOCCERSIZE3D(C_PPDBL4(a,b,c,d))
#define sphere_cap_area_2d(a,b) R_DBL(FUNCNAME_SPHERECAPAREA2D(C_PDBL2(a,b)))
#define sphere_cap_area_3d(a,b) R_DBL(FUNCNAME_SPHERECAPAREA3D(C_PDBL2(a,b)))
#define sphere_cap_area_nd(a,b,c) R_DBL(FUNCNAME_SPHERECAPAREAND(C_PDBL3(a,b,c)))
#define sphere_cap_volume_2d(a,b) R_DBL(FUNCNAME_SPHERECAPVOLUME2D(C_PDBL2(a,b)))
#define sphere_cap_volume_3d(a,b) R_DBL(FUNCNAME_SPHERECAPVOLUME3D(C_PDBL2(a,b)))
#define sphere_distance1(a,b,c,d,e) R_DBL(FUNCNAME_SPHEREDISTANCE1(C_PDBL5(a,b,c,d,e)))
#define sphere_distance2(a,b,c,d,e) R_DBL(FUNCNAME_SPHEREDISTANCE2(C_PDBL5(a,b,c,d,e)))
#define sphere_distance3(a,b,c,d,e) R_DBL(FUNCNAME_SPHEREDISTANCE3(C_PDBL5(a,b,c,d,e)))
#define sphere_exp_contains_point_3d(a,b,c,d,e) R_UCHR(FUNCNAME_SPHEREEXPCONTAINSPOINT3D(C_PPDBL5(a,b,c,d,e)))
#define sphere_exp_point_near_3d(a,b,c,d,e,f) FUNCNAME_SPHEREEXPPOINTNEAR3D(C_PPDBL6(a,b,c,d,e,f))
#define sphere_exp2imp_3d(a,b,c,d,e,f) FUNCNAME_SPHEREEXP2IMP3D(C_PPDBL6(a,b,c,d,e,f))
#define sphere_exp2imp_nd(a,b,c,d) FUNCNAME_SPHEREEXP2IMPND(C_DT3PIT(a,b,c,d))
#define sphere_imp_area_3d(a) R_DBL(FUNCNAME_SPHEREIMPAREA3D(C_SDBL(a)))
#define sphere_imp_area_nd(a,b) R_DBL(FUNCNAME_SPHEREIMPAREAND(C_PDBL2(a,b)))
#define sphere_imp_contains_point_3d(a,b,c) R_UCHR(FUNCNAME_SPHEREIMPCONTAINSPOINT3D(C_IT2PIT(a,b,c)))
#define sphere_imp_grid_icos_size(a,b,c,d) FUNCNAME_SPHEREIMPGRIDICOSSIZE(C_DT3PDT(a,b,c,d))
#define sphere_imp_gridfaces_3d(a,b,c,d,e) FUNCNAME_SPHEREIMPGRIDFACES3D(C_3DT2PDT(a,b,c,d,e))
#define sphere_imp_gridlines_3d(a,b,c,d,e) FUNCNAME_SPHEREIMPGRIDLINES3D(C_3DT2PDT(a,b,c,d,e))
#define sphere_imp_gridpoints_3d(a,b,c,d,e,f,g) FUNCNAME_SPHEREIMPGRIDPOINTS3D(C_ITPIT4DTPDTPIT(a,b,c,d,e,f,g))
#define sphere_imp_gridpoints_icos1(a,b,c) FUNCNAME_SPHEREIMPGRIDPOINTSICOS1(C_2DTPIT(a,b,c))
#define sphere_imp_gridpoints_icos2(a,b,c) FUNCNAME_SPHEREIMPGRIDPOINTSICOS2(C_2DTPIT(a,b,c))
#define sphere_imp_line_project_3d(a,b,c,d,e,f,g,h) R_USHRT(FUNCNAME_SPHEREIMPLINEPROJECT3D(C_ITPITDTPITDTPIT2IT(a,b,c,d,e,f,g,h)))
#define sphere_imp_local2xyz_3d(a,b,c,d,e) FUNCNAME_SPHEREIMPLOCAL2XYZ3D(C_ITPIT2ITPIT(a,b,c,d,e))
#define sphere_imp_point_near_3d(a,b,c,d) FUNCNAME_SPHEREIMPPOINTNEAR3D(C_IT3PIT(a,b,c,d))
#define sphere_imp_point_project_3d(a,b,c,d) FUNCNAME_SPHEREIMPPOINTPROJECT3D(C_IT3PIT(a,b,c,d))
#define sphere_imp_spiralpoints_3d(a,b,c,d) FUNCNAME_SPHEREIMPSPIRALPOINTS3D(C_DTIT2PIT(c,a,b,d))
#define sphere_imp_volume_3d(a) R_DBL(FUNCNAME_SPHEREIMPVOLUME3D(C_SDBL(a)))
#define sphere_imp_volume_nd(a,b) R_DBL(FUNCNAME_SPHEREIMPVOLUMEND(C_DTIT(a,b)))
#define sphere_imp_zone_area_3d(a,b,c) R_DBL(FUNCNAME_SPHEREIMPZONEAREA3D(C_PDBL3(a,b,c)))
#define sphere_imp_zone_volume_3d(a,b,c) R_DBL(FUNCNAME_SPHEREIMPZONEVOLUME3D(C_PDBL3(a,b,c)))
#define sphere_imp2exp_3d(a,b,c,d,e,f) R_DBL(FUNCNAME_SPHEREIMP2EXP3D(C_IT5PIT(a,b,c,d,e,f)))
#define sphere_k(a) R_DBL(FUNCNAME_SPHEREK(C_SUSHRT(a)))
#define sphere_triangle_angles_to_area(a,b,c,d) R_DBL(FUNCNAME_SPHERETRIANGLEANGLESTOAREA(C_PDBL4(a,b,c,d)))
#define sphere_triangle_sides_to_angles(a,b,c,d,e,f,g) FUNCNAME_SPHERETRIANGLESIDESTOANGLES(C_4IT3PIT(a,b,c,d,e,f,g))
#define sphere_triangle_vertices_to_angles(a,b,c,d,e,f,g) FUNCNAME_SPHERETRIANGLEVERTICESTOANGLES(C_IT6PIT(a,b,c,d,e,f,g))
#define sphere_triangle_vertices_to_area(a,b,c,d) R_DBL(FUNCNAME_SPHERETRIANGLEVERTICESTOAREA(C_IT3PIT(a,b,c,d)))
#define sphere_triangle_vertices_to_centroid(a,b,c,d,e) FUNCNAME_SPHERETRIANGLEVERTICESTOCENTROID(C_IT4PIT(a,b,c,d,e))
#define sphere_triangle_vertices_to_orientation(a,b,c) R_SHRT(FUNCNAME_SPHERETRIANGLEVERTICESTOORIENTATION(C_PPDBL3(a,b,c)))
#define sphere_triangle_vertices_to_sides(a,b,c,d,e,f,g) FUNCNAME_SPHERETRIANGLEVERTICESTOSIDES(C_IT6PIT(a,b,c,d,e,f,g))
#define sphere_unit_area_nd(a) R_DBL(FUNCNAME_SPHEREUNITAREAND(C_SUSHRT(a)))
#define sphere_unit_area_values(a,b,c) FUNCNAME_SPHEREUNITAREAVALUES(C_2DTPIT(a,b,c))
#define sphere_unit_sample_2d(a) FUNCNAME_SPHEREUNITSAMPLE2D(a)
#define sphere_unit_sample_3d(a) FUNCNAME_SPHEREUNITSAMPLE3D(a)
#define sphere_unit_sample_3d_2(a) FUNCNAME_SPHEREUNITSAMPLE3D2(a)
#define sphere_unit_sample_nd(a,b) FUNCNAME_SPHEREUNITSAMPLEND(C_DTPI(a,b))
#define sphere_unit_sample_nd_2(a,b) FUNCNAME_SPHEREUNITSAMPLEND2(C_DTPI(a,b))
#define sphere_unit_sample_nd_3(a,b) FUNCNAME_SPHEREUNITSAMPLEND3(C_DTPI(a,b))
#define sphere_unit_volume_nd(a) R_DBL(FUNCNAME_SPHEREUNITVOLUMEND(C_SUSHRT(a)))
#define sphere_unit_volume_values(a,b,c) FUNCNAME_SPHEREUNITVOLUMEVALUES(C_2DTPIT(a,b,c))
#define sphere01_distance_xyz(a,b) R_DBL(FUNCNAME_SPHERE01DISTANCEXYZ(C_PPDBL2(a,b)))
#define sphere01_polygon_area(a,b,c) FUNCNAME_SPHERE01POLYGONAREA(C_DT2PIT(a,b,c))
#define string_2d(a,b,c,d,e,f) R_DBL(FUNCNAME_STRING2D(C_DT2PIT3PDT(a,b,c,d,e,f)))
#define super_ellipse_points_2d(a,b,c,d,e,f,g) FUNCNAME_SUPERELLIPSEPOINT2D(C_4ITDT2PIT(b,c,d,e,f,g,a))
#define tetrahedron_barycentric_3d(a,b) FUNCNAME_TETRAHEDRONBARYCENTRIC3D(C_PPDBL2(a,b))
#define tetrahedron_centroid_3d(a) FUNCNAME_TETRAHEDRONCENTROID3D(a)
#define tetrahedron_circumsphere_3d(a,b,c) R_DBL(FUNCNAME_TETRAHEDRONCIRCUMSPHERE3D(C_PPDBL3(a,b,c)))
#define tetrahedron_contains_point_3d(a,b) R_UCHR(FUNCNAME_TETRAHEDRONCONTAINSPOINT3D(C_PPDBL2(a,b)))
#define tetrahedron_dihedral_angles_3d(a) FUNCNAME_TETRAHEDRONDIHEDRALANGLES3D(a)
#define tetrahedron_edge_length_3d(a) FUNCNAME_TETRAHEDRONEDGELENGTH3D(a)
#define tetrahedron_face_angles_3d(a,b) FUNCNAME_TETRAHEDRONFACEANGLES3D(C_PPDBL2(a,b))
#define tetrahedron_face_areas_3d(a,b) FUNCNAME_TETRAHEDRONFACEAREAS3D(C_PPDBL2(a,b))
#define tetrahedron_insphere_3d(a,b,c) FUNCNAME_TETRAHEDRONINSPHERE3D(C_PPDBL3(a,b,c))
#define tetrahedron_lattice_layer_point_next(a,b,c) FUNCNAME_TETRAHEDRONLATTICELAYERPOINTNEXT(C_PIPDTPB(a,b,c))
#define tetrahedron_lattice_point_next(a,b,c) FUNCNAME_TETRAHEDRONLATTICEPOINTNEXT(C_PIPDTPB(a,b,c))
#define tetrahedron_quality1_3d(a) R_DBL(FUNCNAME_TETRAHEDRONQUALITY13D(a))
#define tetrahedron_quality2_3d(a) R_DBL(FUNCNAME_TETRAHEDRONQUALITY23D(a))
#define tetrahedron_quality3_3d(a) R_DBL(FUNCNAME_TETRAHEDRONQUALITY33D(a))
#define tetrahedron_quality4_3d(a) R_DBL(FUNCNAME_TETRAHEDRONQUALITY43D(a))
#define tetrahedron_rhombic_shape_3d(a,b,c,d,e,f) FUNCNAME_TETRAHEDRONRHOMBICSHAPE3D(C_PIT2DTPDTDTPDT(d,a,b,e,c,f))
#define tetrahedron_rhombic_size_3d(a,b,c,d) FUNCNAME_TETRAHEDRONRHOMBICSIZE3D(C_PPUSHRT4(a,b,c,d))
#define tetrahedron_sample_3d(a,b,c,d) FUNCNAME_TETRAHEDRONSAMPLE3D(C_PITDTPIPIT(a,b,c,d))
#define tetrahedron_shape_3d(a,b,c,d,e,f) FUNCNAME_TETRAHEDRONSHAPE3D(C_3DTPIT2PI(a,b,c,d,e,f))
#define tetrahedron_size_3d(a,b,c,d) FUNCNAME_TETRAHEDRONSIZE3D(C_PPUSHRT4(a,b,c,d))
#define tetrahedron_solid_angles_3d(a) FUNCNAME_TETRAHEDRONSOLIDANGLES3D(a)
#define tetrahedron_unit_lattice_point_num_3d(a) R_USHRT(FUNCNAME_TETRAHEDRONUNITLATTICEPOINTNUM3D(C_SUSHRT(a)))
#define tetrahedron_volume_3d(a) R_DBL(FUNCNAME_TETRAHEDRONVOLUME3D(a))
#define theta2_adjust(a,b) FUNCNAME_THETA2ADJUST(C_PPDBL2(a,b))
#define theta3_adjust(a,b,c) FUNCNAME_THETA3ADJUST(C_PPDBL3(a,b,c))
#define tmat_init(a) FUNCNAME_TMATINIT(a)
#define tmat_mxm(a,b,c) FUNCNAME_TMATMXM(C_PPDBL3(a,b,c))
#define tmat_mxp(a,b,c) FUNCNAME_TMATMXP(C_PPDBL3(a,b,c))
#define tmat_mxp2(a,b,c,d) FUNCNAME_TMATMXP2(C_DT3PIT(b,c,d,a))
#define tmat_mxv(a,b,c) FUNCNAME_TMATMXV(C_PPDBL3(a,b,c))
#define tmat_rot_axis(a,b,c,d) FUNCNAME_TMATROTAXIS(C_2PITITCH(a,b,c,d))
#define tmat_rot_vector(a,b,c,d) FUNCNAME_TMATROTVECTOR(C_IT3PIT(c,a,b,d))
#define tmat_scale(a,b,c) FUNCNAME_TMATSCALE(C_PPDBL3(a,b,c))
#define tmat_shear(a,b,c,d) FUNCNAME_TMATSHEAR(C_2PITPCHIT(a,b,c,d))
#define tmat_trans(a,b,c) FUNCNAME_TMATTRANS(C_PPDBL3(a,b,c))
#define torus_area_3d(a,b) R_DBL(FUNCNAME_TORUSAREA3D(C_PDBL2(a,b)))
#define torus_volume_3d(a,b) R_DBL(FUNCNAME_TORUSVOLUME3D(C_PDBL2(a,b)))
#define tp_to_xyz(a,b) FUNCNAME_TPTOXYZ(C_PDBL2(a,b))
#define triangle_angles_2d(a,b) FUNCNAME_TRIANGLEANGLES2D(C_PPDBL2(a,b))
#define triangle_angles_2d_new(a) FUNCNAME_TRIANGLEANGLES2DNEW(a)
#define triangle_angles_3d(a,b) FUNCNAME_TRIANGLEANGLES3D(C_PPDBL2(a,b))
#define triangle_angles_3d_new(a) FUNCNAME_TRIANGLEANGLES3DNEW(a)
#define triangle_area_2d(a) R_DBL(FUNCNAME_TRIANGLEAREA2D(a))
#define triangle_area_3d(a) R_DBL(FUNCNAME_TRIANGLEAREA3D(a))
#define triangle_area_3d_2(a) R_DBL(FUNCNAME_TRIANGLEAREA3D2(a))
#define triangle_area_3d_3(a) R_DBL(FUNCNAME_TRIANGLEAREA3D3(a))
#define triangle_area_heron(a) R_DBL(FUNCNAME_TRIANGLEAREAHERON(a))
#define triangle_area_vector_3d(a) FUNCNAME_TRIANGLEAREAVECTOR3D(a)
#define triangle_barycentric_2d(a,b) FUNCNAME_TRIANGLEBARYCENTRIC2D(C_PPDBL2(a,b))
#define triangle_centroid_2d(a) FUNCNAME_TRIANGLECENTROID2D(a)
#define triangle_centroid_3d(a) FUNCNAME_TRIANGLECENTROID3D(a)
#define triangle_circumcenter_2d(a) FUNCNAME_TRIANGLECIRCUMCENTER2D(a)
#define triangle_circumcenter_2d_2(a) FUNCNAME_TRIANGLECIRCUMCENTER2D2(a)
#define triangle_circumcenter(a,b) FUNCNAME_TRIANGLECIRCUMCENTER(C_DTPIT(a,b))
#define triangle_circumcircle_2d(a,b,c) FUNCNAME_TRIANGLECIRCUMCIRCLE2D(C_PPDBL3(a,b,c))
#define triangle_circumcircle_2d_2(a,b,c) FUNCNAME_TRIANGLECIRCUMCIRCLE2D2(C_PPDBL3(a,b,c))
#define triangle_circumradius_2d(a) R_DBL(FUNCNAME_TRIANGLECIRCUMRADIUS2D(a))
#define triangle_contains_line_exp_3d(a,b,c,d,e) FUNCNAME_TRIANGLECONTAINSLINEEXP3D(C_3PITPBPIT(a,b,c,d,e))
#define triangle_contains_line_par_3d(a,b,c,d,e) FUNCNAME_TRIANGLECONTAINSLINEPAR3D(C_3PITPBPIT(a,b,c,d,e))
#define triangle_contains_point_2d_1(a,b) R_UCHR(FUNCNAME_TRIANGLECONTAINSPOINT2D1(C_PPDBL2(a,b)))
#define triangle_contains_point_2d_2(a,b) R_UCHR(FUNCNAME_TRIANGLECONTAINSPOINT2D2(C_PPDBL2(a,b)))
#define triangle_contains_point_2d_3(a,b) R_UCHR(FUNCNAME_TRIANGLECONTAINSPOINT2D3(C_PPDBL2(a,b)))
#define triangle_diameter_2d(a) R_DBL(FUNCNAME_TRIANGLEDIAMETER2D(a))
#define triangle_edge_length_2d(a) FUNCNAME_TRIANGLEEDGELENGTH2D(a)
#define triangle_gridpoints_2d(a,b,c,d,e) FUNCNAME_TRIANGLEGRIDPOINTS2D(C_PIT2DTPDTPIT(a,b,c,d,e))
#define triangle_incenter_2d(a,b) FUNCNAME_TRIANGLEINCENTER2D(C_PPDBL2(a,b))
#define triangle_incircle_2d(a,b,c) FUNCNAME_TRIANGLEINCIRCLE2D(C_PPDBL3(a,b,c))
#define triangle_inradius_2d(a) R_DBL(FUNCNAME_TRIANGLEINRADIUS2D(a))
#define triangle_is_degenerate_nd(a,b) R_UCHR(FUNCNAME_TRIANGLEISDEGENERATEND(C_DTPIT(a,b)))
#define triangle_lattice_layer_point_next(a,b,c) FUNCNAME_TRIANGLELATTICELAYERPOINTNEXT(C_PIPDTPB(a,b,c))
#define triangle_lattice_point_next(a,b,c) FUNCNAME_TRIANGLELATTICEPOINTNEXT(C_PIPDTPB(a,b,c))
#define triangle_line_imp_int_2d(a,b,c,d,e) FUNCNAME_TRIANGLELINEIMPINT2D(C_PIT3ITPDTPIT(a,b,c,d,e))
#define triangle_orientation_2d(a) R_USHRT(FUNCNAME_TRIANGLEORIENTATION2D(a))
#define triangle_orthocenter_2d(a,b,c) FUNCNAME_TRIANGLEORTHOCENTER2D(C_2PITPB(a,b,c))
#define triangle_point_dist_2d(a,b) R_DBL(FUNCNAME_TRIANGLEPOINTDIST2D(C_PPDBL2(a,b)))
#define triangle_point_dist_3d(a,b) R_DBL(FUNCNAME_TRIANGLEPOINTDIST3D(C_PPDBL2(a,b)))
#define triangle_point_dist_signed_2d(a,b) R_DBL(FUNCNAME_TRIANGLEPOINTDISTSIGNED2D(C_PPDBL2(a,b)))
#define triangle_point_near_2d(a,b,c,d) FUNCNAME_TRIANGLEPOINTNEAR2D(C_PPDBL4(a,b,c,d))
#define triangle_quality_2d(a) R_DBL(FUNCNAME_TRIANGLEQUALITY2D(a))
#define triangle_right_lattice_point_num_2d(a,b) R_USHRT(FUNCNAME_TRIANGLERIGHTLATTICEPOINTNUM2D(C_PUSHRT2(a,b)))
#define triangle_sample(a,b,c,d) FUNCNAME_TRIANGLESAMPLE(C_PITDTPIPIT(a,b,c,d))
#define triangle_unit_lattice_point_num_2d(a) R_USHRT(FUNCNAME_TRIANGLEUNITLATTICEPOINTNUM2D(C_SUSHRT(a)))
#define triangle_xsi_to_xy_2d(a,b,c) FUNCNAME_TRIANGLEXSITOXY2D(C_PPDBL3(a,b,c))
#define triangle_xy_to_xsi_2d(a,b,c) FUNCNAME_TRIANGLEXYTOXSI2D(C_PPDBL3(a,b,c))
#define truncated_octahedron_shape_3d(a,b,c,d,e,f) FUNCNAME_TRUNCATEDOCTAHEDRONSHAPE3D(C_3DTPIT2PI(a,b,c,d,e,f))
#define truncated_octahedron_size_3d(a,b,c,d) FUNCNAME_TRUNCATEDOCTAHEDRONSIZE3D(C_PPUSHRT4(a,b,c,d))
#define tube_2d(a,b,c,d,e) FUNCNAME_TUBE2D(C_DTIT3PIT(b,a,c,d,e))
#define tuple_next2(a,b,c,d,e) FUNCNAME_TUPLENEXT2(C_DT4PDT(a,b,c,d,e))
#define vector_directions_nd(a,b,c) FUNCNAME_VECTORDIRECTIONSND(C_DT2PIT(a,b,c))
#define vector_rotate_2d(a,b,c) FUNCNAME_VECTORROTATE2D(C_IT2PIT(b,a,c))
#define vector_rotate_3d(a,b,c,d) FUNCNAME_VECTORROTATE3D(C_IT3PIT(c,a,b,d))
#define vector_rotate_base_2d(a,b,c,d) FUNCNAME_VECTORROTATEBASE2D(C_IT3PIT(c,a,b,d))
#define vector_separation_2d(a,b) R_DBL(FUNCNAME_VECTORSEPARATION2D(C_PPDBL2(a,b)))
#define vector_separation_3d(a,b) R_DBL(FUNCNAME_VECTORSEPARATION3D(C_PPDBL2(a,b)))
#define vector_separation_nd(a,b,c) R_DBL(FUNCNAME_VECTORSEPARATIONND(C_DT2PIT(a,b,c)))
#define vector_unit_nd(a,b) FUNCNAME_VECTORUNITND(C_DTPIT(a,b))
#define voxels_dist_l1_3d(a,b) R_USHRT(FUNCNAME_VOXELSDISTL13D(C_PPUSHRT2(a,b)))
#define voxels_dist_l1_nd(a,b,c) R_USHRT(FUNCNAME_VOXELSDISTL1ND(C_DT2PIT(a,b,c)))
#define voxels_line_3d(a,b,c,d) FUNCNAME_VOXELSLINE3D(C_DT3PDT(c,a,b,d))
#define voxels_region_3d(a,b,c,d,e,f,g,h,i) FUNCNAME_VOXELSREGION3D(C_4DT4PDT(a,b,c,d,e,f,g,h,i))
#define voxels_step_3d(a,b,c,d,e,f) FUNCNAME_VOXELSSTEP3D(C_PDTPI3DTPI(a,b,c,d,e,f))
#define xy_to_polar(a,b,c) FUNCNAME_XYTOPOLAR(C_PPDBL3(a,b,c))
#define xyz_to_radec(a,b,c) FUNCNAME_XYZTORADEC(C_PPDBL3(a,b,c))
#define xyz_to_rtp(a,b,c,d) FUNCNAME_XYZTORTP(C_PPDBL4(a,b,c,d))
#define xyz_to_tp(a,b,c) FUNCNAME_XYZTOTP(C_PPDBL3(a,b,c))

__MATHSUITE __JBURKARDT void * FUNCNAME_DQRANK(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DQRLS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DQRLSS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATAMAX(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATSOLVE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ANGLEBOX2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ANGLECONTAINSRAY2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ANGLERAD2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ANGLERAD3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ANGLERADND(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ANGLETURN2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ANGLEDEG2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ANGLEIDEG2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ANGLEIRAD2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BOX01CONTAINSPOINT2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BOX01CONTAINSPOINTND(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BOXCONTAINSPOINT2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BOXCONTAINSPOINTND(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BOXRAYINT2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BOXSEGMENTCLIP2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CIRCLEARCPOINTNEAR2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CIRCLEDIA2IMP2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CIRCLEDIA2IMP3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CIRCLEEXPCONTAINSPOINT2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CIRCLEEXP2IMP2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CIRCLEIMPCONTAINSPOINT2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CIRCLEIMPLINEPARINT2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CIRCLEIMPPOINTDIST2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CIRCLEIMPPOINTDISTSIGNED2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CIRCLEIMPPOINTNEAR2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CIRCLEIMPPOINTSARC2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CIRCLEIMP2EXP2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CIRCLELUNEAREA2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CIRCLEPPPR2IMP3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CIRCLESECTORAREA2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CIRCLESECTORCONTAINSPOINT2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CIRCLETRIANGLEAREA2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CIRCLETRIPLEANGLES2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CIRCLESIMPINT2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CONV3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CUBESIZE3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CUBESHAPE3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CYLINDERPOINTDIST3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CYLINDERPOINTDISTSIGNED3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CYLINDERPOINTINSIDE3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CYLINDERVOLUME3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DGEDET(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DGE0FA(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DGE0SL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DEGREESTORADIANS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DMSTORADIANS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DISKPOINTDIST3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DODECSHAPE3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DODECSIZE3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DUALSHAPE3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DUALSIZE3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ELLIPSEPOINTDIST2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ELLIPSEPOINTS2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ELLIPSEPOINTSARC2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ELLIPSEAREA2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ENORM0ND(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GLOB2LOC3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HAVERSINE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HALFPLANECONTAINSPOINT2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HALFSPACEIMPTTRIANGLEINT3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HALFSPACENORMTRIANGLEINT3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HALFSPACETRIANGLEINT3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HELIXSHAPE3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HEXAGONCONTAINSPOINT2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HEXAGONSHAPE2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HEXAGONVERTICES2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_I4FACTORIAL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_I4FACTORIAL2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_I4SWAP(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ICOSSHAPE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ICOSSIZE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LINEEXPISDEGENERATEND(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LINEEXPPOINTDIST2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LINEEXPPOINTDIST3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LINEEXPPOINTDISTSIGNED2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LINEEXPPOINTNEAR2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LINEEXPPOINTNEAR3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LINEEXP2IMP2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LINEEXP2PAR2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LINEEXP2PAR3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LINEIMPPOINTDIST2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LINEIMPPOINTDISTSIGNED2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LINEIMP2EXP2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LINEIMP2PAR2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LINEPARPOINTDIST2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LINEPARPOINTDIST3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LINEPAR2EXP2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LINEPAR2EXP3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LINEPAR2IMP2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LINESEXPANGLE3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LINESEXPANGLEND(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LINESEXPDIST3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LINESEXPDIST3D2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LINESEXPEQUAL2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LINESEXPINT2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LINESEXPNEAR3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LINESEXPPARALLEL2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LINESEXPPARALLEL3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LINESIMPINT2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LINESPARINT2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LOC2GLOB3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_MINABS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_MINQUAD(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_OCTAHEDRONSHAPE3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_OCTAHEDRONSIZE3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PARABOLAEX(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PARABOLAEX2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PARALLELOGRAMAREA2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PARALLELOGRAMAREA3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PARALLELOGRAMCONTAINSPOINT2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PARALLELOGRAMCONTAINSPOINT3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PARALLELOGRAMPOINTDIST3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PARALLEPIPEDCONTAINSPOINT3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PARALLEPIPEDPOINTDIST3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PLANEEXPGRID3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PLANEEXPPOINTDIST3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PLANEEXPNORMAL3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PLANEEXPPRO2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PLANEEXPPRO3(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PLANEEXPPROJECT3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PLANEEXP2IMP3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PLANEEXP2NORMAL3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PLANEIMPLINEPARINT3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PLANEIMPPOINTDIST3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PLANEIMPPOINTDISTSIGNED3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PLANEIMPPOINTNEAR3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PLANEIMPSEGMENTNEAR3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PLANEIMPTRIANGLEINT3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PLANEIMPTRIANGLEINTADD3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PLANEIMPTRIANGLENEAR3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PLANEIMP2EXP3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PLANEIMP2NORMAL3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PLANENORMALBASIS3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PLANENORMALLINEEXPINT3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PLANENORMALTETRAHEDRONINTERSECT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PLANENORMALTRIANGLEINT3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PLANENORMALUNIFORM3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PLANENORMALUNIFORMND(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PLANENORMAL2EXP3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PLANENORMAL2IMP3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PLANESIMPANGLE3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POINTSAVOIDPOINTNAIVE2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POINTSBISECTLINEIMP2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POINTSBISECTLINEPAR2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POINTSCENTROID2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POINTSCOLIN2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POINTSCOLIN3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POINTSDIST2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POINTSDIST3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POINTSDISTND(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POINTSHULL2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POINTSPOINTNEARNAIVE2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POINTSPOINTNEARNAIVE3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POINTSPOINTNEARNAIVEND(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POLARTOXY(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POLYGON12D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POLYGONAREA2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POLYGONAREA2D2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POLYGONAREA3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POLYGONAREA3D2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POLYGONCONTAINSPOINT2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POLYGONCONTAINSPOINT2D2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POLYGONDIAMETER2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POLYGONINRADDATA2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POLYGONISCONVEX(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POLYGONOUTRADDATA2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POLYGONSIDEDATA2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POLYGONSOLIDANGLE3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POLYGONX2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POLYGONY2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POLYGONXX2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POLYGONXY2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POLYGONYY2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POLYGONLATTICEAREA2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POLYHEDRONAREA3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POLYHEDRONCONTAINSPOINT3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POLYHEDRONVOLUME3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POLYHEDRONVOLUME3D2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POLYLINELENGTHND(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PROVEC(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PYRAMIDVOLUME3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_QUADAREA2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_QUADAREA2D2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_QUADAREA3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_QUADCONTAINSPOINT2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_QUADPOINTDIST2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_QUADPOINTDISTSIGNED2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_QUADPOINTNEAR2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_QUATNORM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECANGLE3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECBRACKET(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_RADIANSTODMS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_RADIANSTODEGREES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ROTATIONAXISVECTOR3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ROTATIONAXIS2MAT3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ROTATIONAXIS2QUAT3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ROTATIONMATVECTOR3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ROTATIONMAT2AXIS3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ROTATIONMAT2QUAT3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ROTATIONQUATVECTOR3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ROTATIONQUAT2AXIS3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ROTATIONQUAT2MAT3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_RTPTOXYZ(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SEGMENTCONTAINSPOINT1D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SEGMENTCONTAINSPOINT2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SEGMENTPOINTCOORDS2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SEGMENTPOINTCOORDS3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SEGMENTPOINTDIST2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SEGMENTPOINTDIST3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SEGMENTPOINTNEAR2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SEGMENTPOINTNEAR3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SEGMENTSCURVATURE2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SEGMENTSDIST2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SEGMENTSDIST3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SEGMENTSDIST3DOLD(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SEGMENTSINT1D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SEGMENTSINT2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SHAPEPOINTDIST2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SHAPEPOINTNEAR2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SHAPERAYINT2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SIMPLEXLATTICELAYERPOINTNEXT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SIMPLEXLATTICEPOINTNEXT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SIMPLEXVOLUMEND(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SOCCERSHAPE3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SOCCERSIZE3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHEREDIA2IMP3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHEREEXPCONTAINSPOINT3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHEREEXPPOINTNEAR3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHEREEXP2IMP3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHEREEXP2IMPND(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHEREIMPCONTAINSPOINT3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHEREIMPGRIDICOSSIZE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHEREIMPGRIDFACES3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHEREIMPGRIDLINES3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHEREIMPGRIDPOINTS3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHEREIMPGRIDPOINTSICOS1(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHEREIMPGRIDPOINTSICOS2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHEREIMPLINEPROJECT3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHEREIMPLOCAL2XYZ3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHEREIMPPOINTNEAR3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHEREIMPPOINTPROJECT3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHEREIMPSPIRALPOINTS3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHEREIMP2EXP3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHERETRIANGLEANGLESTOAREA(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHERETRIANGLESIDESTOANGLES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHERETRIANGLEVERTICESTOANGLES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHERETRIANGLEVERTICESTOAREA(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHERETRIANGLEVERTICESTOCENTROID(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHERETRIANGLEVERTICESTOORIENTATION(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHERETRIANGLEVERTICESTOSIDES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHERE01DISTANCEXYZ(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHERE01POLYGONAREA(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_STRING2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SUPERELLIPSEPOINT2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TETRAHEDRONCIRCUMSPHERE3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TETRAHEDRONCONTAINSPOINT3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TETRAHEDRONFACEANGLES3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TETRAHEDRONFACEAREAS3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TETRAHEDRONINSPHERE3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TETRAHEDRONLATTICELAYERPOINTNEXT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TETRAHEDRONLATTICEPOINTNEXT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TETRAHEDRONQUALITY13D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TETRAHEDRONQUALITY23D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TETRAHEDRONQUALITY33D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TETRAHEDRONQUALITY43D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TETRAHEDRONRHOMBICSHAPE3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TETRAHEDRONRHOMBICSIZE3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TETRAHEDRONSAMPLE3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TETRAHEDRONSHAPE3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TETRAHEDRONSIZE3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TETRAHEDRONVOLUME3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TETRAHEDRONUNITLATTICEPOINTNUM3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_THETA2ADJUST(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_THETA3ADJUST(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TMATINIT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TMATMXM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TMATMXP(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TMATMXP2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TMATMXV(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TMATROTAXIS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TMATROTVECTOR(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TMATSCALE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TMATSHEAR(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TMATTRANS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLEANGLES2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLEANGLES3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLEAREA2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLEAREA3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLEAREA3D2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLEAREA3D3(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLEAREAHERON(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLECIRCUMCIRCLE2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLECIRCUMCIRCLE2D2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLECIRCUMRADIUS2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLECONTAINSLINEEXP3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLECONTAINSLINEPAR3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLECONTAINSPOINT2D1(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLECONTAINSPOINT2D2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLECONTAINSPOINT2D3(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLEDIAMETER2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLEGRIDPOINTS2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLEINCENTER2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLEINCIRCLE2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLEINRADIUS2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLEISDEGENERATEND(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLELATTICELAYERPOINTNEXT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLELATTICEPOINTNEXT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLELINEIMPINT2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLEORIENTATION2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLEORTHOCENTER2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLEPOINTDIST2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLEPOINTDIST3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLEPOINTDISTSIGNED2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLEPOINTNEAR2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLEQUALITY2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLERIGHTLATTICEPOINTNUM2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLESAMPLE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLEXSITOXY2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLEXYTOXSI2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRUNCATEDOCTAHEDRONSHAPE3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRUNCATEDOCTAHEDRONSIZE3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TUBE2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TUPLENEXT2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_VECTORDIRECTIONSND(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_VECTORROTATE2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_VECTORROTATE3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_VECTORROTATEBASE2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_VECTORSEPARATION2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_VECTORSEPARATION3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_VECTORSEPARATIONND(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_VECTORUNITND(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_VOXELSDISTL13D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_VOXELSDISTL1ND(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_VOXELSLINE3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_VOXELSREGION3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_VOXELSSTEP3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_XYTOPOLAR(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_XYZTORADEC(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_XYZTORTP(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_XYZTOTP(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_QRSOLVE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ANGLEHALF2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ANNULUSSECTORCENTROID2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BALLUNITSAMPLE2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BALLUNITSAMPLE3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BALLUNITSAMPLEND(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BASISMAP3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CIRCLEIMPPOINTS2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CIRCLEIMPPOINTS3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CIRCLELLR2IMP2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CIRCLELUNECENTROID2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CIRCLEPPR2IMP2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CIRCLESECTORCENTROID2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CONECENTROID3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CYLINDERPOINTNEAR3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CYLINDERSAMPLE3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DIRECTIONPERT3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DIRECTIONUNIFORM2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DIRECTIONUNIFORM3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DIRECTIONUNIFORMND(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ELLIPSEPOINTNEAR2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LINEEXPNORMAL2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LINEEXPPERP2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LINEPARPOINTNEAR2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LINEPARPOINTNEAR3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PLANENORMALQRTOXYZ(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PLANENORMALXYZTOQR(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POINTSPOINTSNEARNAIVE2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POINTSPOINTSNEARNAIVE3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POLYGONANGLES2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POLYGONCENTROID2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POLYGONCENTROID2D2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POLYGONCENTROID3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POLYGONEXPAND2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POLYGONNORMAL3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POLYHEDRONCENTROID3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POLYLINEARCLENGTHND(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POLYLINEINDEXPOINTND(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POLYLINEPOINTSND(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POLYLOOPARCLENGTHND(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POLYLOOPPOINTSND(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_QUATCONJ(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_QUATINV(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_QUATMUL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATMVNEW(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8MATSOLVE2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_RADECTOXYZ(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_RADECDISTANCE3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHEREUNITSAMPLE2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHEREUNITSAMPLE3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHEREUNITSAMPLE3D2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHEREUNITSAMPLEND(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHEREUNITSAMPLEND2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHEREUNITSAMPLEND3(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHEREUNITVOLUMEND(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHEREUNITVOLUMEVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TETRAHEDRONBARYCENTRIC3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TETRAHEDRONCENTROID3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TETRAHEDRONDIHEDRALANGLES3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TETRAHEDRONEDGELENGTH3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TETRAHEDRONSOLIDANGLES3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TPTOXYZ(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLEANGLES2DNEW(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLEANGLES3DNEW(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLEAREAVECTOR3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLEBARYCENTRIC2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLECENTROID2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLECENTROID3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLECIRCUMCENTER2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLECIRCUMCENTER2D2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLECIRCUMCENTER(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLEEDGELENGTH2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GETSEED(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HEXAGONUNITAREA2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LINESPARANGLE3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LINESPARDIST3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CIRCLEAREA2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_COTRAD(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_HEXAGONAREA2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SIMPLEXUNITVOLUMEND(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHEREIMPAREA3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHEREIMPVOLUME3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHEREK(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHEREUNITAREAND(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHEREUNITVOLUMEND(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHEREUNITAREAVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLEUNITLATTICEPOINTNUM2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ANNULUSAREA2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ANNULUSSECTORAREA2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CONEAREA3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CONEVOLUME3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_I4DEDEKINDFACTOR(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_I4DEDEKINDSUM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_I4GCD(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_I4MODP(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_I4SIGN(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SIMPLEXUNITLATTICEPOINTND(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHERECAPAREA2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHERECAPAREA3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHERECAPVOLUME2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHERECAPVOLUME3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHEREIMPAREAND(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHEREIMPVOLUMEND(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TORUSAREA3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TORUSVOLUME3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LINEIMPISDEGENERATE2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PLANEIMPISDEGENERATE3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHERECAPAREAND(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHERECAPVOLUMEND(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHEREIMPZONEAREA3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHEREIMPZONEVOLUME3D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHEREDISTANCE1(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHEREDISTANCE2(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHEREDISTANCE3(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LINESIMPANGLE2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LINESIMPDIST2D(void *);

#endif
