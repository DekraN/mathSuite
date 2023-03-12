#ifndef __DISABLEDEEP_SPHEREQUAD

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  ityp   sphere01_quad_icos1c ( const register dim_typ factor,void fun ( dim_typ n, ityp x[], ityp v[] ), int *node_num )
/******************************************************************************/
/*
  Purpose:
    SPHERE01_QUAD_ICOS1C: centroid rule, subdivide then project.
  Discussion:
    This function estimates an integral over the surface of the unit sphere.
    This function sets up an icosahedral grid, and subdivides each
    edge of the icosahedron into FACTOR subedges.  These edges define a grid
    within each triangular icosahedral face.  The centroids of these
    triangles can be determined.  All of these calculations are done,
    essentially, on the FLAT faces of the icosahedron.  Only then are
    the triangle vertices and centroids projected to the sphere.
    The resulting grid of spherical triangles and projected centroids
    is used to apply a centroid quadrature rule over the surface of
    the unit sphere.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 September 2010
  Author:
    John Burkardt
  Parameters:
    Input, int FACTOR, the subdivision factor, which must
    be at least 1.
    Input, void FUN ( int n, double x[], double v[] ), evaluates the
    integrand.
    Output, int *NODE_NUM, the number of evaluation points.
    Output, double SPHERE01_QUAD_ICOS1C, the estimated integral.
*/
{
    dim_typ a;
    ityp a_xyz[3];
    ityp *a2_xyz;
    ityp area;
    ityp area_total;
    dim_typ b;
    ityp b_xyz[3];
    ityp *b2_xyz;
    dim_typ c;
    ityp c_xyz[3];
    ityp *c2_xyz;
    dim_typ edge_num;
    int *edge_point;
    dim_typ f1;
    dim_typ f2;
    dim_typ f3;
    dim_typ face;
    dim_typ face_num;
    int *face_order;
    int *face_point;
    dim_typ face_order_max;
    dim_typ i;
    ityp *node_xyz;
    ityp *point_coord;
    dim_typ point_num;
    ityp result;
    ityp v[1];
    /*
    Size the icosahedron.
    */
    icos_size ( &point_num, &edge_num, &face_num, &face_order_max );
    /*
    Set the icosahedron.
    */
    point_coord = ( ityp * ) malloc ( 3 * point_num * sizeof ( ityp ) );
    edge_point = ( int * ) malloc ( edge_num * sizeof ( int ) <<1 );
    face_order = ( int * ) malloc ( face_num * sizeof ( int ) );
    face_point = ( int * ) malloc ( face_order_max * face_num * sizeof ( int ) );

    icos_shape ( point_num, edge_num, face_num, face_order_max,
    point_coord, edge_point, face_order, face_point );
    /*
    Initialize the integral data.
    */
    result = area_total = 0.00;
    *node_num = 0;
    /*
    Pick a face of the icosahedron, and identify its vertices as A, B, C.
    */
    for ( face = 0; face < face_num; ++face )
    {
        a = face_point[0+face*3];
        b = face_point[1+face*3];
        c = face_point[2+face*3];

        #pragma omp parallel for num_threads(3)
        for ( i = 0; i < 3; ++i )
        {
            a_xyz[i] = point_coord[i+a*3];
            b_xyz[i] = point_coord[i+b*3];
            c_xyz[i] = point_coord[i+c*3];
        }
        /*
        Some subtriangles will have the same direction as the face.
        Generate each in turn, by determining the barycentric coordinates
        of the centroid (F1,F2,F3), from which we can also work out the barycentric
        coordinates of the vertices of the subtriangle.
        */
        for ( f3 = 1; f3 <= 3 * factor - 2; f3 += 3 )
        {
            for ( f2 = 1; f2 <= 3 * factor - f3 - 1; f2 += 3 )
            {
                f1 = 3 * factor - f3 - f2;

                node_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1, f2, f3 );

                a2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1 + 2, f2 - 1, f3 - 1 );
                b2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1 - 1, f2 + 2, f3 - 1 );
                c2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1 - 1, f2 - 1, f3 + 2 );

                area = sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz );

                fun ( 1, node_xyz, v );

                ++ *node_num;
                result += area * v[0];
                area_total += area;

                free ( a2_xyz );
                free ( b2_xyz );
                free ( c2_xyz );
                free ( node_xyz );
            }
        }
        /*
        The other subtriangles have the opposite direction from the face.
        Generate each in turn, by determining the barycentric coordinates
        of the centroid (F1,F2,F3), from which we can also work out the barycentric
        coordinates of the vertices of the subtriangle.
        */
        for ( f3 = 2; f3 <= 3 * factor - 4; f3 += 3 )
        {
            for ( f2 = 2; f2 <= 3 * factor - f3 - 2; f2 += 3 )
            {
                f1 = 3 * factor - f3 - f2;

                node_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1, f2, f3 );

                a2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1 - 2, f2 + 1, f3 + 1 );
                b2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1 + 1, f2 - 2, f3 + 1 );
                c2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1 + 1, f2 + 1, f3 - 2 );

                area = sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz );

                fun ( 1, node_xyz, v );

                ++ *node_num;
                result += area * v[0];
                area_total += area;

                free ( a2_xyz );
                free ( b2_xyz );
                free ( c2_xyz );
                free ( node_xyz );
            }
        }
    }
    /*
    Discard allocated memory.
    */
    free ( edge_point );
    free ( face_order );
    free ( face_point );
    free ( point_coord );

    return result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  ityp   sphere01_quad_icos1m ( const register dim_typ factor,void fun ( dim_typ n, ityp x[], ityp v[] ), int *node_num )
/******************************************************************************/
/*
  Purpose:
    SPHERE01_QUAD_ICOS1M: midside rule, subdivide then project.
  Discussion:
    This function estimates an integral over the surface of the unit sphere.
    This function sets up an icosahedral grid, and subdivides each
    edge of the icosahedron into FACTOR subedges.  These edges define a grid
    within each triangular icosahedral face.  The midsides of these
    triangles can be determined.  All of these calculations are done,
    essentially, on the FLAT faces of the icosahedron.  Only then are
    the triangle vertices and midsides projected to the sphere.
    The resulting grid of spherical triangles and projected midsides
    is used to apply a midside quadrature rule over the surface of
    the unit sphere.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    24 September 2010
  Author:
    John Burkardt
  Parameters:
    Input, int FACTOR, the subdivision factor, which must
    be at least 1.
    Input, void FUN ( int n, double x[], double v[] ), evaluates the
    integrand.
    Output, int *NODE_NUM, the number of evaluation points.
    Output, double SPHERE01_QUAD_ICOS1M, the estimated integral.
*/
{
    dim_typ a;
    ityp a_xyz[3];
    ityp *a2_xyz;
    ityp *a3_xyz;
    ityp area;
    ityp area_total;
    dim_typ b;
    ityp b_xyz[3];
    ityp *b2_xyz;
    ityp *b3_xyz;
    dim_typ c;
    ityp c_xyz[3];
    ityp *c2_xyz;
    ityp *c3_xyz;
    dim_typ edge;
    dim_typ edge_num;
    int *edge_point;
    dim_typ f;
    dim_typ f1;
    dim_typ f2;
    dim_typ f3;
    dim_typ face;
    dim_typ face_num;
    int *face_order;
    int *face_point;
    dim_typ face_order_max;
    dim_typ i;
    dim_typ j;
    dim_typ node;
    ityp node_norm;
    ityp *node_xyz;
    ityp *point_coord;
    dim_typ point_num;
    ityp result;
    ityp va[1];
    ityp vb[1];
    ityp vc[1];
    /*
    Size the icosahedron.
    */
    icos_size ( &point_num, &edge_num, &face_num, &face_order_max );
    /*
    Set the icosahedron.
    */
    point_coord = ( ityp * ) malloc ( 3 * point_num * sizeof ( ityp ) );
    edge_point = ( int * ) malloc ( edge_num * sizeof ( int ) <<1 );
    face_order = ( int * ) malloc ( face_num * sizeof ( int ) );
    face_point = ( int * ) malloc ( face_order_max * face_num * sizeof ( int ) );

    icos_shape ( point_num, edge_num, face_num, face_order_max,
    point_coord, edge_point, face_order, face_point );
    /*
    Initialize the integral data.
    */
    result = area_total = 0.00;
    *node_num = 0;
    /*
    Pick a face of the icosahedron, and identify its vertices as A, B, C.
    */
    for ( face = 0; face < face_num; ++face )
    {
    a = face_point[0+face*3];
        b = face_point[1+face*3];
        c = face_point[2+face*3];

        #pragma omp parallel for num_threads(3)
        for ( i = 0; i < 3; ++i )
        {
            a_xyz[i] = point_coord[i+a*3];
            b_xyz[i] = point_coord[i+b*3];
            c_xyz[i] = point_coord[i+c*3];
        }
        /*
        Deal with subtriangles that have same orientation as face.
        */
        for ( f1 = 0; f1 <= factor - 1; ++f1)
        {
            for ( f2 = 0; f2 <= factor - f1 - 1; ++f2 )
            {
                f3 = factor - f1 - f2;

                a2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1 + 1, f2,     f3 - 1 );
                b2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1,     f2 + 1, f3 - 1 );
                c2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1,     f2,     f3 );

                area = sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz );

                a3_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, (f1<<1) + 1, (f2<<1) + 1, (f3<<1) - 2 );
                b3_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, (f1<<1),  (f2<<1) + 1, (f3<<1) - 1 );
                c3_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, (f1<<1) + 1, (f2<<1),  (f3<<1) - 1 );

                *node_num += 3;
                fun ( 1, a3_xyz, va );
                fun ( 1, b3_xyz, vb );
                fun ( 1, c3_xyz, vc );
                result += area * ( va[0] + vb[0] + vc[0] ) / 3.00;
                area_total += area;

                free ( a2_xyz );
                free ( a3_xyz );
                free ( b2_xyz );
                free ( b3_xyz );
                free ( c2_xyz );
                free ( c3_xyz );
            }
        }
        /*
        Deal with subtriangles that have opposite orientation as face.
        */
        for ( f3 = 0; f3 <= factor - 2; ++f3 )
        {
            for ( f2 = 1; f2 <= factor - f3 - 1; ++f2)
            {
                f1 = factor - f2 - f3;

                a2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1 - 1, f2,     f3 + 1 );
                b2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1,     f2 - 1, f3 + 1 );
                c2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1,     f2,     f3 );

                area = sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz );

                a3_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, (f1<<1) - 1, (f2<<1) - 1, (f3<<1) + 2 );
                b3_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, (f1<<1),  (f2<<1) - 1, (f3<<1) + 1 );
                c3_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, (f1<<1) - 1, (f2<<1),  (f3<<1) + 1 );

                *node_num += 3;
                fun ( 1, a3_xyz, va );
                fun ( 1, b3_xyz, vb );
                fun ( 1, c3_xyz, vc );
                result += area * ( va[0] + vb[0] + vc[0] ) / 3.00;
                area_total += area;

                free ( a2_xyz );
                free ( a3_xyz );
                free ( b2_xyz );
                free ( b3_xyz );
                free ( c2_xyz );
                free ( c3_xyz );
            }
        }
    }
    /*
    Discard allocated memory.
    */
    free ( edge_point );
    free ( face_order );
    free ( face_point );
    free ( point_coord );

    return result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  ityp   sphere01_quad_icos1v ( const register dim_typ factor,void fun ( dim_typ n, ityp x[], ityp v[] ), int *node_num )
/******************************************************************************/
/*
  Purpose:
    SPHERE01_QUAD_ICOS1V: vertex rule, subdivide then project.
  Discussion:
    This function estimates an integral over the surface of the unit sphere.
    This function sets up an icosahedral grid, and subdivides each
    edge of the icosahedron into FACTOR subedges.  These edges define a grid
    within each triangular icosahedral face.  The vertices of these
    triangles can be determined.  All of these calculations are done,
    essentially, on the FLAT faces of the icosahedron.  Only then are
    the triangle vertices projected to the sphere.
    The resulting grid of spherical triangles is used to apply a vertex
    quadrature rule over the surface of the unit sphere.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 September 2010
  Author:
    John Burkardt
  Parameters:
    Input, int FACTOR, the subdivision factor, which must
    be at least 1.
    Input, void FUN ( int n, double x[], double v[] ), evaluates the
    integrand.
    Output, int *NODE_NUM, the number of evaluation points.
    Output, double SPHERE01_QUAD_ICOS2M, the estimated integral.
*/
{
    dim_typ a;
    ityp a_xyz[3];
    ityp *a2_xyz;
    ityp area;
    ityp area_total;
    dim_typ b;
    ityp b_xyz[3];
    ityp *b2_xyz;
    dim_typ c;
    ityp c_xyz[3];
    ityp *c2_xyz;
    dim_typ edge;
    dim_typ edge_num;
    int *edge_point;
    dim_typ f;
    dim_typ f1;
    dim_typ f2;
    dim_typ f3;
    dim_typ face;
    dim_typ face_num;
    int *face_order;
    int *face_point;
    dim_typ face_order_max;
    dim_typ i;
    dim_typ j;
    dim_typ node;
    ityp *point_coord;
    dim_typ point_num;
    ityp result;
    ityp va[1];
    ityp vb[1];
    ityp vc[1];
    /*
    Size the icosahedron.
    */
    icos_size ( &point_num, &edge_num, &face_num, &face_order_max );
    /*
    Set the icosahedron.
    */
    point_coord = ( ityp * ) malloc ( 3 * point_num * sizeof ( ityp ) );
    edge_point = ( int * ) malloc ( edge_num * sizeof ( int ) <<1 );
    face_order = ( int * ) malloc ( face_num * sizeof ( int ) );
    face_point = ( int * ) malloc ( face_order_max * face_num * sizeof ( int ) );

    icos_shape ( point_num, edge_num, face_num, face_order_max,
    point_coord, edge_point, face_order, face_point );
    /*
    Initialize the integral data.
    */
    result = area_total = 0.00;
    *node_num = 0;

    /*
    Pick a face of the icosahedron, and identify its vertices as A, B, C.
    */
    for ( face = 0; face < face_num; ++face )
    {
        a = face_point[0+face*3];
        b = face_point[1+face*3];
        c = face_point[2+face*3];

        #pragma omp parallel for num_threads(3)
        for ( i = 0; i < 3; ++i )
        {
            a_xyz[i] = point_coord[i+a*3];
            b_xyz[i] = point_coord[i+b*3];
            c_xyz[i] = point_coord[i+c*3];
        }
        /*
        Deal with subtriangles that have same orientation as face.
        */
        for ( f1 = 0; f1 <= factor - 1; ++f1 )
        {
            for ( f2 = 0; f2 <= factor - f1 - 1; ++f2 )
            {
                f3 = factor - f1 - f2;

                a2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1 + 1, f2,     f3 - 1 );
                b2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1,     f2 + 1, f3 - 1 );
                c2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1,     f2,     f3 );

                area = sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz );

                *node_num += 3;
                fun ( 1, a2_xyz, va );
                fun ( 1, b2_xyz, vb );
                fun ( 1, c2_xyz, vc );
                result += area * ( va[0] + vb[0] + vc[0] ) / 3.00;
                area_total += area;

                free ( a2_xyz );
                free ( b2_xyz );
                free ( c2_xyz );
            }
        }
        /*
        Deal with subtriangles that have opposite orientation as face.
        */
        for ( f3 = 0; f3 <= factor - 2; ++f3)
        {
            for ( f2 = 1; f2 <= factor - f3 - 1; ++f2 )
            {
                f1 = factor - f2 - f3;

                a2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1 - 1, f2,     f3 + 1 );
                b2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1,     f2 - 1, f3 + 1 );
                c2_xyz = sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1,     f2,     f3 );

                area = sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz );

                *node_num += 3;
                fun ( 1, a2_xyz, va );
                fun ( 1, b2_xyz, vb );
                fun ( 1, c2_xyz, vc );
                result += area * ( va[0] + vb[0] + vc[0] ) / 3.00;
                area_total += area;

                free ( a2_xyz );
                free ( b2_xyz );
                free ( c2_xyz );
            }
        }
    }
    /*
    Discard allocated memory.
    */
    free ( edge_point );
    free ( face_order );
    free ( face_point );
    free ( point_coord );

    return result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  ityp   sphere01_quad_icos2v ( const register dim_typ factor,void fun ( dim_typ n, ityp x[], ityp v[] ), int *node_num )
/******************************************************************************/
/*
  Purpose:
    SPHERE01_QUAD_ICOS2V: vertex rule, subdivide then project.
  Discussion:
    This function estimates an integral over the surface of the unit sphere.
    This function sets up an icosahedral grid, and subdivides each
    edge of the icosahedron into FACTOR subedges.  These edges define a grid
    within each triangular icosahedral face.  The vertices of these
    triangles can be determined.  All of these calculations are done,
    essentially, on the FLAT faces of the icosahedron.  Only then are
    the triangle vertices projected to the sphere.
    The resulting grid of spherical triangles is used to apply a vertex
    quadrature rule over the surface of the unit sphere.
    This is a revision of SPHERE01_QUAD_ICOS2V that attempted to use a more
    sophisticated scheme to map points from the planar triangle to the surface
    of the unit sphere.  Very little improvement to the estimated integral
    was observed, so development of this scheme has been set aside for now.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 September 2010
  Author:
    John Burkardt
  Parameters:
    Input, int FACTOR, the subdivision factor, which must
    be at least 1.
    Input, void FUN ( int n, double x[], double v[] ), evaluates the
    integrand.
    Output, int *NODE_NUM, the number of evaluation points.
    Output, double SPHERE01_QUAD_ICOS2V, the estimated integral.
*/
{
    dim_typ a;
    ityp a_xyz[3];
    ityp *a2_xyz;
    ityp area;
    ityp area_total;
    dim_typ b;
    ityp b_xyz[3];
    ityp *b2_xyz;
    dim_typ c;
    ityp c_xyz[3];
    ityp *c2_xyz;
    dim_typ edge;
    dim_typ edge_num;
    int *edge_point;
    dim_typ f;
    dim_typ f1;
    dim_typ f2;
    dim_typ f3;
    dim_typ face;
    dim_typ face_num;
    int *face_order;
    int *face_point;
    dim_typ face_order_max;
    dim_typ i;
    dim_typ j;
    dim_typ node;
    ityp *point_coord;
    dim_typ point_num;
    ityp result;
    ityp va[1];
    ityp vb[1];
    ityp vc[1];
    /*
    Size the icosahedron.
    */
    icos_size ( &point_num, &edge_num, &face_num, &face_order_max );
    /*
    Set the icosahedron.
    */
    point_coord = ( ityp * ) malloc ( 3 * point_num * sizeof ( ityp ) );
    edge_point = ( int * ) malloc ( edge_num * sizeof ( int ) <<1 );
    face_order = ( int * ) malloc ( face_num * sizeof ( int ) );
    face_point = ( int * ) malloc ( face_order_max * face_num * sizeof ( int ) );

    icos_shape ( point_num, edge_num, face_num, face_order_max,
    point_coord, edge_point, face_order, face_point );
    /*
    Initialize the integral data.
    */
    result = area_total = 0.00;
    *node_num = 0;
    /*
    Pick a face of the icosahedron, and identify its vertices as A, B, C.
    */
    for ( face = 0; face < face_num; ++face)
    {
        a = face_point[0+face*3];
        b = face_point[1+face*3];
        c = face_point[2+face*3];

        #pragma omp parallel for num_threads(3)
        for ( i = 0; i < 3; ++i )
        {
            a_xyz[i] = point_coord[i+a*3];
            b_xyz[i] = point_coord[i+b*3];
            c_xyz[i] = point_coord[i+c*3];
        }
        /*
        Deal with subtriangles that have same orientation as face.
        */
        for ( f1 = 0; f1 <= factor - 1; ++f1 )
        {
            for ( f2 = 0; f2 <= factor - f1 - 1; ++f2)
            {
            f3 = factor - f1 - f2;

            a2_xyz = sphere01_triangle_project2 ( a_xyz, b_xyz, c_xyz, f1 + 1, f2,     f3 - 1 );
            b2_xyz = sphere01_triangle_project2 ( a_xyz, b_xyz, c_xyz, f1,     f2 + 1, f3 - 1 );
            c2_xyz = sphere01_triangle_project2 ( a_xyz, b_xyz, c_xyz, f1,     f2,     f3 );

            area = sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz );

            *node_num += 3;
            fun ( 1, a2_xyz, va );
            fun ( 1, b2_xyz, vb );
            fun ( 1, c2_xyz, vc );
            result = result + area * ( va[0] + vb[0] + vc[0] ) / 3.00;
            area_total += area;

            free ( a2_xyz );
            free ( b2_xyz );
            free ( c2_xyz );
            }
        }
        /*
        Deal with subtriangles that have opposite orientation as face.
        */
        for ( f3 = 0; f3 <= factor - 2; ++f3 )
        {
            for ( f2 = 1; f2 <= factor - f3 - 1; ++f2 )
            {
                f1 = factor - f2 - f3;

                a2_xyz = sphere01_triangle_project2 ( a_xyz, b_xyz, c_xyz, f1 - 1, f2,     f3 + 1 );
                b2_xyz = sphere01_triangle_project2 ( a_xyz, b_xyz, c_xyz, f1,     f2 - 1, f3 + 1 );
                c2_xyz = sphere01_triangle_project2 ( a_xyz, b_xyz, c_xyz, f1,     f2,     f3 );

                area = sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz );

                *node_num += 3;
                fun ( 1, a2_xyz, va );
                fun ( 1, b2_xyz, vb );
                fun ( 1, c2_xyz, vc );
                result += area * ( va[0] + vb[0] + vc[0] ) / 3.00;
                area_total += area;

                free ( a2_xyz );
                free ( b2_xyz );
                free ( c2_xyz );
            }
        }
    }
    /*
    Discard allocated memory.
    */
    free ( edge_point );
    free ( face_order );
    free ( face_point );
    free ( point_coord );

    return result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  ityp   sphere01_quad_llc ( void f ( dim_typ n, ityp x[], ityp v[] ), const register ityp h,dim_typ *n )
/******************************************************************************/
/*
  Purpose:
    SPHERE01_QUAD_LLC: Longitude/Latitude grid with centroid rule.
  Discussion:
    The sphere is broken up into spherical triangles, whose sides
    do not exceed the length H.  Then a centroid rule is used on
    each spherical triangle.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 September 2010
  Author:
    John Burkardt.
  Parameters:
    Input, void F ( int n, double x[], double v[] ), evaluates the
    integrand.
    Input, double H, the maximum length of a side of the spherical
    quadrilaterals.
    Output, int *N, the number of points used.
    Output, double SPHERE01_QUAD_LLC, the approximate integral.
*/
{
    ityp area;
    dim_typ i, j;
    ityp phi;
    int phi_num;
    ityp phi1;
    ityp phi2;
    ityp result;
    ityp sector_area;
    ityp sphere_area;
    ityp theta;
    dim_typ theta_num;
    ityp theta1;
    ityp theta2;
    ityp v[1];
    ityp *x;
    ityp *x1;
    ityp *x11;
    ityp *x12;
    ityp *x2;
    ityp *x21;
    ityp *x22;
    /*
    Choose PHI and THETA counts that make short sides.
    */
    phi_num = ( int ) ( M_PI / h );

    if ( h * ( ityp ) ( phi_num ) < M_PI )
        ++ phi_num;

    theta_num = ( int ) ( M_2TPI / h );

    if ( h * ( ityp ) ( theta_num ) < M_PI )
        ++ theta_num;

    *n = 0;
    result = 0.00;
    /*
    Only one THETA (and hence, only one PHI.)
    */
    if ( theta_num == 1 )
    {
        sphere_area = 4.00 * M_PI;

        theta = 0.00;
        phi = M_PI / 2.00;
        x = tp_to_xyz ( theta, phi );

        f ( 1, x, v );
        ++ *n;
        result = sphere_area * v[0];
        free ( x );
    }
    /*
    Several THETA''s, one PHI.
    */
    else if ( phi_num == 1 )
    {
        sphere_area = 4.00 * M_PI;
        sector_area = sphere_area / ( ityp ) ( theta_num );

        result = 0.00;

        for ( j = 1; j <= theta_num; ++j )
        {
            theta = ( ityp ) ( ( j - 1 ) <<1 ) * M_PI / ( ityp ) ( theta_num );
            phi = M_PI / 2.00;
            x = tp_to_xyz ( theta, phi );
            f ( 1, x, v );
            ++ *n;
            result += sector_area * v[0];
            free ( x );
        }
    }
    /*
    At least two PHI''s.
    */
    else
    {
        result = 0.00;
        /*
        Picture in top row, with V1 = north pole:

        V1
        /  \
        /    \
        V12----V22
        */
        phi1 = 0.00;
        phi2 = M_PI / ( ityp ) ( phi_num );

        for ( j = 1; j <= theta_num; ++j )
        {
            theta1 = ( ityp ) ( j - 1 ) * M_2TPI / ( ityp ) ( theta_num );
            theta2 = ( ityp ) ( j ) * M_2TPI / ( ityp ) ( theta_num );

            x1 = tp_to_xyz ( theta1, phi1 );
            x12 = tp_to_xyz ( theta1, phi2 );
            x22 = tp_to_xyz ( theta2, phi2 );

            area = sphere01_triangle_vertices_to_area ( x1, x12, x22 );
            x = sphere01_triangle_vertices_to_centroid ( x1, x12, x22 );
            f ( 1, x, v );
            ++ *n;
            result += area * v[0];
            free ( x );
            free ( x1 );
            free ( x12 );
            free ( x22 );
        }
        /*
        Picture in all intermediate rows:

        V11--V21
        | \  |
        |  \ |
        V12--V22
        */
        for ( i = 2; i <= phi_num - 1; ++i)
        {
            phi1 = ( ityp ) ( i - 1 ) * M_PI / ( ityp ) ( phi_num );
            phi2 = ( ityp ) ( i ) * M_PI / ( ityp ) ( phi_num );

            for ( j = 1; j <= theta_num; ++j )
            {
                theta1 = ( ityp ) ( j - 1 ) * M_2TPI / ( ityp ) ( theta_num );
                theta2 = ( ityp ) ( j ) * M_2TPI / ( ityp ) ( theta_num );

                x11 = tp_to_xyz ( theta1, phi1 );
                x21 = tp_to_xyz ( theta2, phi1 );
                x12 = tp_to_xyz ( theta1, phi2 );
                x22 = tp_to_xyz ( theta2, phi2 );

                area = sphere01_triangle_vertices_to_area ( x11, x12, x22 );
                x = sphere01_triangle_vertices_to_centroid ( x11, x12, x22 );
                f ( 1, x, v );
                ++ *n;
                result +=  area * v[0];
                free ( x );

                area = sphere01_triangle_vertices_to_area ( x22, x21, x11 );
                x = sphere01_triangle_vertices_to_centroid ( x22, x21, x11 );
                f ( 1, x, v );
                ++ *n;
                result += area * v[0];
                free ( x );
                free ( x11 );
                free ( x12 );
                free ( x21 );
                free ( x22 );
            }
        }
        /*
        Picture in last row, with V2 = south pole:

        V11----V21
        \    /
        \  /
        V2
        */
        phi1 = ( ityp ) ( phi_num - 1 ) * M_PI / ( ityp ) ( phi_num );
        phi2 = M_PI;

        for ( j = 1; j <= theta_num; ++j )
        {
            theta1 = ( ityp ) ( j - 1 ) * M_2TPI / ( ityp ) ( theta_num );
            theta2 = ( ityp ) ( j ) * M_2TPI / ( ityp ) ( theta_num );

            x11 = tp_to_xyz ( theta1, phi1 );
            x21 = tp_to_xyz ( theta2, phi1 );
            x2 = tp_to_xyz ( theta2, phi2 );

            area = sphere01_triangle_vertices_to_area ( x11, x2, x21 );
            x = sphere01_triangle_vertices_to_centroid ( x11, x2, x21 );
            f ( 1, x, v );
            ++ *n;
            result += area * v[0];
            free ( x );
            free ( x11 );
            free ( x2 );
            free ( x21 );
        }
    }
    return result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  ityp   sphere01_quad_llm ( void f ( dim_typ n, ityp x[], ityp v[] ), const register ityp h,dim_typ *n )
/******************************************************************************/
/*
  Purpose:
    SPHERE01_QUAD_LLM: longitude/latitude grid plus midside rule.
  Discussion:
    The sphere is broken up into spherical triangles, whose sides
    do not exceed the length H.  Then the function is evaluated
    at the midsides, and the average is multiplied by the area.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 September 2010
  Author:
    John Burkardt
  Parameters:
    Input, void F ( int n, double x[], double v[] ), evaluates the
    integrand.
    Input, double H, the maximum length of a side of the spherical
    quadrilaterals.
    Output, int *N, the number of points used.
    Output, double SPHERE01_QUAD_LLM, the approximate integral.
*/
{
    ityp area;
    dim_typ i, j;
    ityp m1[3];
    ityp m2[3];
    ityp m3[3];
    ityp phi;
    dim_typ phi_num;
    ityp phi1;
    ityp phi2;
    ityp result;
    ityp sector_area;
    ityp sphere_area;
    ityp theta;
    dim_typ theta_num;
    ityp theta1;
    ityp theta2;
    ityp v[1];
    ityp *x;
    ityp *x1;
    ityp *x11;
    ityp *x12;
    ityp *x2;
    ityp *x21;
    ityp *x22;
    /*
    Choose PHI and THETA counts that make short sides.
    */
    phi_num = ( int ) ( M_PI / h );

    if ( h * ( ityp ) ( phi_num ) < M_PI )
        ++ phi_num;

    theta_num = ( int ) ( M_2TPI / h );

    if ( h * ( ityp ) ( theta_num ) < M_PI )
        ++ theta_num;

    *n = 0;
    result = 0.00;
    /*
    Only one THETA (and hence, only one PHI.)
    */
    if ( theta_num == 1 )
    {
        sphere_area = 4.00 * M_PI;

        theta = 0.00;
        phi = M_PI / 2.00;
        x = tp_to_xyz ( theta, phi );
        f ( 1, x, v );
        ++ *n;
        result = sphere_area * v[0];
        free ( x );
    }
    /*
    Several THETA''s, one PHI.
    */
    else if ( phi_num == 1 )
    {
        sphere_area = 4.00 * M_PI;
        sector_area /= ( ityp ) ( theta_num );

        result = 0.0;

        for ( j = 1; j <= theta_num; ++j )
        {
            theta = ( ityp ) ( ( j - 1 ) <<1 ) * M_PI / ( ityp ) ( theta_num );
            phi = M_PI / 2.00;
            x = tp_to_xyz ( theta, phi );
            f ( 1, x, v );
            ++ *n;
            result += sector_area * v[0];
            free ( x );
        }
    }
    /*
    At least two PHI''s.
    */
    else
    {
        result = phi1 = 0.00;
        /*
        Picture:

        V1
        /  \
        /    \
        V12----V22
        */
        phi2 = M_PI / ( ityp ) ( phi_num );

        for ( j = 1; j <= theta_num; ++j )
        {
            theta1 = ( ityp ) ( j - 1 ) * M_2TPI / ( ityp ) ( theta_num );
            theta2 = ( ityp ) ( j ) * M_2TPI / ( ityp ) ( theta_num );

            x1 = tp_to_xyz ( theta1, phi1 );
            x12 = tp_to_xyz ( theta1, phi2 );
            x22 = tp_to_xyz ( theta2, phi2 );

            area = sphere01_triangle_vertices_to_area ( x1, x12, x22 );

            sphere01_triangle_vertices_to_midpoints ( x1, x12, x22, m1, m2, m3 );

            f ( 1, m1, v );
            ++ *n;
            result += area * v[0] / 3.00;
            f ( 1, m2, v );
            ++ *n;
            result += area * v[0] / 3.00;
            f ( 1, m3, v );
            ++ *n;
            result += area * v[0] / 3.00;
            free ( x1 );
            free ( x12 );
            free ( x22 );
        }
        /*
        Picture:

        V11--V21
        | \  |
        |  \ |
        V12--V22
        */
        for ( i = 2; i <= phi_num - 1; ++i )
        {
            phi1 = ( ityp ) ( i - 1 ) * M_PI / ( ityp ) ( phi_num );
            phi2 = ( ityp ) ( i ) * M_PI / ( ityp ) ( phi_num );

            for ( j = 1; j <= theta_num; ++j )
            {
                theta1 = ( ityp ) ( j - 1 ) * M_2TPI / ( ityp ) ( theta_num );
                theta2 = ( ityp ) ( j ) * M_2TPI / ( ityp ) ( theta_num );

                x11 = tp_to_xyz ( theta1, phi1 );
                x21 = tp_to_xyz ( theta2, phi1 );
                x12 = tp_to_xyz ( theta1, phi2 );
                x22 = tp_to_xyz ( theta2, phi2 );

                area = sphere01_triangle_vertices_to_area ( x11, x12, x22 );

                sphere01_triangle_vertices_to_midpoints ( x11, x12, x22, m1, m2, m3 );

                f ( 1, m1, v );
                ++ *n;
                result += area * v[0] / 3.00;
                f ( 1, m2, v );
                ++ *n;
                result += area * v[0] / 3.00;
                f ( 1, m3, v );
                ++ *n;
                result += area * v[0] / 3.00;

                area = sphere01_triangle_vertices_to_area ( x22, x21, x11 );

                sphere01_triangle_vertices_to_midpoints ( x22, x21, x11, m1, m2, m3 );

                f ( 1, m1, v );
                ++ *n;
                result += area * v[0] / 3.00;
                f ( 1, m2, v );
                ++ *n;
                result += area * v[0] / 3.0;
                f ( 1, m3, v );
                ++ *n;
                result += area * v[0] / 3.00;
                free ( x11 );
                free ( x12 );
                free ( x21 );
                free ( x22 );
            }
        }
        /*
        Picture:

        V11----V21
        \    /
        \  /
        V2
        */
        phi1 = ( ityp ) ( phi_num - 1 ) * M_PI / ( ityp ) ( phi_num );
        phi2 = M_PI;

        for ( j = 1; j <= theta_num; ++j )
        {
            theta1 = ( ityp ) ( j - 1 ) * M_2TPI / ( ityp ) ( theta_num );
            theta2 = ( ityp ) ( j ) * M_2TPI / ( ityp ) ( theta_num );

            x11 = tp_to_xyz ( theta1, phi1 );
            x21 = tp_to_xyz ( theta2, phi1 );
            x2 = tp_to_xyz ( theta2, phi2 );

            area = sphere01_triangle_vertices_to_area ( x11, x2, x21 );

            sphere01_triangle_vertices_to_midpoints ( x11, x2, x21, m1, m2, m3 );

            f ( 1, m1, v );
            ++ *n;
            result += area * v[0] / 3.00;
            f ( 1, m2, v );
            ++ *n;
            result += area * v[0] / 3.00;
            f ( 1, m3, v );
            ++ *n;
            result += area * v[0] / 3.00;
            free ( x11 );
            free ( x2 );
            free ( x21 );
        }
    }

    return result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  ityp   sphere01_quad_llv ( void f ( dim_typ n, ityp x[], ityp v[] ), const register ityp h,dim_typ *n )
/******************************************************************************/
/*
  Purpose:
    SPHERE01_QUAD_LLV: longitude/latitude grid with vertex rule.
  Discussion:
    The sphere is broken up into spherical triangles, whose sides
    do not exceed the length H.  Then the function is evaluated
    at the vertices, and the average is multiplied by the area.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 September 2010
  Author:
    John Burkardt
  Parameters:
    Input, void F ( int n, double x[], double v[] ), evaluates the
    integrand.
    Input, double H, the maximum length of a side of the spherical
    quadrilaterals.
    Output, int *N, the number of points used.
    Output, double SPHERE01_QUAD_LLV, the approximate integral.
*/
{
    ityp area;
    dim_typ i, j;
    ityp phi;
    dim_typ phi_num;
    ityp phi1;
    ityp phi2;

    ityp result;
    ityp sector_area;
    ityp sphere_area;
    ityp theta;
    dim_typ theta_num;
    ityp theta1;
    ityp theta2;
    ityp v[1];
    ityp *x;
    ityp *x1;
    ityp *x11;
    ityp *x12;
    ityp *x2;
    ityp *x21;
    ityp *x22;
    /*
    Choose PHI and THETA counts that make short sides.
    */
    phi_num = ( int ) ( M_PI / h );

    if ( h * ( ityp ) ( phi_num ) < M_PI )
        ++ phi_num;

    theta_num = ( int ) ( M_2TPI / h );

    if ( h * ( ityp ) ( theta_num ) < M_PI )
        ++ theta_num;

    *n = 0;
    result = 0.00;
    /*
    Only one THETA (and hence, only one PHI.)
    */
    if ( theta_num == 1 )
    {
        sphere_area = 4.00 * M_PI;

        theta = 0.00;
        phi = M_PI / 2.00;
        x = tp_to_xyz ( theta, phi );
        f ( 1, x, v );
        result = sphere_area * v[0];
        free ( x );
    }
    /*
    Several THETA''s, one PHI.
    */
    else if ( phi_num == 1 )
    {
        sphere_area = 4.00 * M_PI;
        sector_area /= ( ityp ) ( theta_num );

        result = 0.00;

        for ( j = 1; j <= theta_num; ++j )
        {
            theta = ( ityp ) ( ( j - 1 ) <<1 ) * M_PI / ( ityp ) ( theta_num );
            phi = M_PI / 2.00;
            x = tp_to_xyz ( theta, phi );
            f ( 1, x, v );
            *n = *n + 1;
            result += sector_area * v[0];
            free ( x );
        }
    }
    /*
    At least two PHI''s.
    */
    else
    {
        result = 0.00;
        /*
        Picture:

        V1
        /  \
        /    \
        V12----V22
        */
        phi1 = 0.00;
        phi2 = M_PI / ( ityp ) ( phi_num );

        for ( j = 1; j <= theta_num; ++j )
        {
            theta1 = ( ityp ) ( j - 1 ) * M_2TPI / ( ityp ) ( theta_num );
            theta2 = ( ityp ) ( j ) * M_2TPI / ( ityp ) ( theta_num );

            x1 = tp_to_xyz ( theta1, phi1 );
            x12 = tp_to_xyz ( theta1, phi2 );
            x22 = tp_to_xyz ( theta2, phi2 );

            area = sphere01_triangle_vertices_to_area ( x1, x12, x22 );

            f ( 1, x1, v );
            *n = *n + 1;
            result += area * v[0] / 3.00;
            f ( 1, x12, v );
            *n = *n + 1;
            result += area * v[0] / 3.00;
            f ( 1, x22, v );
            *n = *n + 1;
            result += area * v[0] / 3.00;
            free ( x1 );
            free ( x12 );
            free ( x22 );
        }
        /*
        Picture:

        V11--V21
        | \  |
        |  \ |
        V12--V22
        */
        for ( i = 2; i <= phi_num - 1; ++i )
        {
            phi1 = ( ityp ) ( i - 1 ) * M_PI / ( ityp ) ( phi_num );
            phi2 = ( ityp ) ( i ) * M_PI / ( ityp ) ( phi_num );

            for ( j = 1; j <= theta_num; ++j )
            {
                theta1 = ( ityp ) ( j - 1 ) * M_2TPI / ( ityp ) ( theta_num );
                theta2 = ( ityp ) ( j ) * M_2TPI / ( ityp ) ( theta_num );

                x11 = tp_to_xyz ( theta1, phi1 );
                x21 = tp_to_xyz ( theta2, phi1 );
                x12 = tp_to_xyz ( theta1, phi2 );
                x22 = tp_to_xyz ( theta2, phi2 );

                area = sphere01_triangle_vertices_to_area ( x11, x12, x22 );

                f ( 1, x11, v );
                ++ *n;
                result += area * v[0] / 3.00;
                f ( 1, x12, v );
                ++ *n;
                result += area * v[0] / 3.00;
                f ( 1, x22, v );
                ++ *n;
                result += area * v[0] / 3.00;

                area = sphere01_triangle_vertices_to_area ( x22, x21, x11 );

                f ( 1, x22, v );
                ++ *n;
                result += area * v[0] / 3.00;
                f ( 1, x21, v );
                ++ *n;
                result += area * v[0] / 3.00;
                f ( 1, x11, v );
                ++ *n;
                result += area * v[0] / 3.00;

                free ( x11 );
                free ( x12 );
                free ( x21 );
                free ( x22 );
            }
        }
        /*
        Picture:

        V11----V21
        \    /
        \  /
        V2
        */
        phi1 = ( ityp ) ( phi_num - 1 ) * M_PI / ( ityp ) ( phi_num );
        phi2 = M_PI;

        for ( j = 1; j <= theta_num; j++ )
        {
            theta1 = ( ityp ) ( j - 1 ) * M_2TPI / ( ityp ) ( theta_num );
            theta2 = ( ityp ) ( j ) * M_2TPI / ( ityp ) ( theta_num );

            x11 = tp_to_xyz ( theta1, phi1 );
            x21 = tp_to_xyz ( theta2, phi1 );
            x2 = tp_to_xyz ( theta2, phi2 );

            area = sphere01_triangle_vertices_to_area ( x11, x2, x21 );

            f ( 1, x11, v );
            ++ *n;
            result += area * v[0] / 3.00;
            f ( 1, x2, v );
            ++ *n;
            result += area * v[0] / 3.00;
            f ( 1, x21, v );
            ++ *n;
            result += area * v[0] / 3.00;

            free ( x11 );
            free ( x2 );
            free ( x21 );
        }
    }
    return result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  ityp   sphere01_quad_mc ( void f ( dim_typ n, ityp x[], ityp v[] ), const register ityp h,int *seed, const register dim_typ n )
/******************************************************************************/
/*
  Purpose:
    SPHERE01_QUAD_MC uses the Monte Carlo rule for sphere quadrature.
  Discussion:
    A number of points N are chosen at random on the sphere, with N
    being determined so that, if the points were laid out on a regular
    grid, the average spacing would be no more than H.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 September 2010
  Author:
    John Burkardt
  Parameters:
    Input, void F ( int n, double x[], double v[] ), evaluates the
    integrand.
    Input, double H, the maximum length of a side of the spherical
    quadrilaterals.
    Input/output, int *SEED, a seed for the random number generator.
    Input, int N, the number of points used.
    Output, double SPHERE01_QUAD_MC, the approximate integral.
*/
{
    ityp result;
    ityp sphere_area;
    ityp *v;
    ityp *x;
    sphere_area = 4.00 * M_PI;
    x = sphere01_sample ( n, seed );
    v = ( ityp * ) malloc ( n * sizeof ( ityp ) );
    f ( n, x, v );
    result = sphere_area * r8vec_sum ( n, v ) / ( ityp ) ( n );
    free ( v );
    free ( x );
    return result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _sphere01_quad_mc_size ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE01_QUAD_MC_SIZE sizes a Monte Carlo rule for sphere quadrature.
  Discussion:
    A number of points N are chosen at random on the sphere, with N
    being determined so that, if the points were laid out on a regular
    grid, the average spacing would be no more than H.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 September 2010
  Author:
    John Burkardt
  Parameters:
    Input, double H, the maximum length of a side of the spherical
    quadrilaterals.
    Output, int SPHERE01_QUAD_MC_SIZE, the number of points to use.
*/
{
	static dim_typ result = USHRT_MAX;
	
	const register ityp h = *(ityp *) data;
	
	result = MAX ( ( int ) ( 4.00 * M_PI / h / h ), 1 );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _sphere01_triangle_angles_to_area ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE01_TRIANGLE_ANGLES_TO_AREA: area of a spherical triangle on the unit sphere.
  Discussion:
    A unit sphere centered at 0 in 3D satisfies the equation:
      X^2 + Y^2 + Z^2 = 1
    A spherical triangle is specified by three points on the surface
    of the sphere.
    The area formula is known as Girard's formula.
    The area of a spherical triangle on the unit sphere is:
      AREA = A + B + C - M_PI
    where A, B and C are the (surface) angles of the triangle.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 September 2010
  Author:
    John Burkardt
  Parameters:
    Input, double A, B, C, the angles of the triangle.
    Output, double SPHERE01_TRIANGLE_ANGLES_TO_AREA, the area of the spherical triangle.
*/
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	const register ityp a = a_data[0];
	const register ityp b = a_data[1];
	const register ityp c = a_data[2];
	
	result = a + b + c - M_PI; 
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _sphere01_triangle_project ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE01_TRIANGLE_PROJECT projects from a plane triangle to a spherical triangle.
  Discussion:
    We assume that points A, B and C lie on the unit sphere, and they
    thus define a spherical triangle.
    They also, of course, define a planar triangle.
    Let (F1,F2,F3) be the barycentric coordinates of a point in this
    planar triangle.
    This function determines the coordinates of the point in the planar
    triangle identified by the barycentric coordinates, and returns the
    coordinates of the projection of that point onto the unit sphere.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 September 2010
  Author:
    John Burkardt
  Parameters:
    Input, double A_XYZ[3], B_XYZ[3], C_XYZ[3], the coordinates
    of the points A, B, and C.
    Input, int F1, F2, F3, the barycentric coordinates
    of a point in the triangle ABC.  Normally, these coordinates would
    be real numbers, and would sum to 1.  For convenience, we allow these
    to be integers which must be divided by F1+F2+F3.
    Output, double NODE_XYZ[3], the coordinates of the
    point on the unit sphere which is the projection of the point on the plane
    whose barycentric coordinates with respect to A, B, and C is
 (F1,F2,F3)/(F1+F2+F3).
*/
{
	const _3dt3pit * const s_data = data;
	
	const register dim_typ f1 = s_data->a0;
	const register dim_typ f2 = s_data->a1;
	const register dim_typ f3 = s_data->a2;
	ityp * a_xyz = s_data->a3;
	ityp * b_xyz = s_data->a4;
	ityp * c_xyz = s_data->a5;
	
	
    dim_typ i;
    ityp *node_xyz;
    ityp norm;

    node_xyz = ( ityp * ) malloc ( 3 * sizeof ( ityp ) );

    #pragma omp parallel for num_threads(3)
    for ( i = 0; i < 3; ++i )
        node_xyz[i] =( ( ityp ) ( f1           ) * a_xyz[i]+ ( ityp ) (      f2      ) * b_xyz[i]+ ( ityp ) (           f3 ) * c_xyz[i] )/ ( ityp ) ( f1 + f2 + f3 );
    norm = r8vec_norm ( 3, node_xyz );

    #pragma omp parallel for num_threads(3)
    for ( i = 0; i < 3; ++i)
        node_xyz[i] /= norm;
    return node_xyz;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _sphere01_triangle_project2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE01_TRIANGLE_PROJECT2 projects from a plane triangle to a spherical triangle.
  Discussion:
    We assume that points A, B and C lie on the unit sphere, and they
    thus define a spherical triangle.
    They also, of course, define a planar triangle.
    Let (F1,F2,F3) be the barycentric coordinates of a point in this
    planar triangle.
    This function determines the coordinates of the point in the planar
    triangle identified by the barycentric coordinates, and returns the
    coordinates of the projection of that point onto the unit sphere.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 September 2010
  Author:
    John Burkardt
  Parameters:
    Input, double A_XYZ(3), B_XYZ(3), C_XYZ(3), the coordinates
    of the points A, B, and C.
    Input, int F1, F2, F3, the barycentric coordinates
    of a point in the triangle ABC.  Normally, these coordinates would
    be real numbers, and would sum to 1.  For convenience, we allow these
    to be integers which must be divided by F1+F2+F3.
    Output, double SPHERE01_TRIANGLE_PROJECT2[3], the coordinates of the
    point on the unit sphere which is the projection of the point on the
    plane whose barycentric coordinates with respect to A, B, and C is
 (F1,F2,F3)/(F1+F2+F3).
*/
{
	const _3dt3pit * const s_data = data;
	
	const register dim_typ f1 = s_data->a0;
	const register dim_typ f2 = s_data->a1;
	const register dim_typ f3 = s_data->a2;
	ityp * a_xyz = s_data->a3;
	ityp * b_xyz = s_data->a4;
	ityp * c_xyz = s_data->a5;
	
    ityp ab[3];
    ityp ac[3];
    ityp acn[3];
    ityp acp[3];
    ityp angle;
    ityp bn[3];
    ityp bp[3];
    ityp cn[3];
    ityp cp[3];
    dim_typ i;
    ityp *node_xyz;
    ityp norm;
    ityp theta_ab;
    ityp theta_ac;
    ityp theta_bc;

    node_xyz = ( ityp * ) malloc ( 3 * sizeof ( ityp ) );
    /*
    This check avoids 0/0 calculations later.
    */
    if ( f2 == 0 && f3 == 0 )
    {
        #pragma omp parallel for num_threads(3)
        for ( i = 0; i < 3; ++i )
            node_xyz[i] = a_xyz[i];
        return node_xyz;
    }
    else if ( f1 == 0 && f3 == 0 )
    {
        #pragma omp parallel for num_threads(3)
        for ( i = 0; i < 3; ++i )
            node_xyz[i] = b_xyz[i];
        return node_xyz;
    }
    else if ( f1 == 0 && f2 == 0 )
    {
        #pragma omp parallel for num_threads(3)
        for ( i = 0; i < 3; ++i )
            node_xyz[i] = c_xyz[i];
        return node_xyz;
    }
    /*
    Determine the angular distances (A,B) and (A,C).
    */
    theta_ab = sphere01_distance_xyz ( a_xyz, b_xyz );
    theta_ac = sphere01_distance_xyz ( a_xyz, c_xyz );
    /*
    Polarize B = BP + BN
    Normalize BN,
    Same for C.
    */
    r8vec_polarize ( 3, b_xyz, a_xyz, bn, bp );
    norm = r8vec_norm ( 3, bn );
    for ( i = 0; i < 3; ++i )
    bn[i] /= norm;
    r8vec_polarize ( 3, c_xyz, a_xyz, cn, cp );
    norm = r8vec_norm ( 3, cn );
    for ( i = 0; i < 3; ++i)
    cn[i] /= norm;
    /*
    Determine AB and AC that use cos ( ( F2 + F3 ) / ( F1 + F2 + F3 ) ) of A
    and cos ( F1 / ( F1 + F2 + F3 ) ) of B or C.
    */
    angle = ( ( ityp ) ( f2 + f3 ) * theta_ab ) / ( ityp ) ( f1 + f2 + f3 );
    #pragma omp parallel for num_threads(3)
    for ( i = 0; i < 3; ++i )
        ab[i] = cos ( angle ) * a_xyz[i] + sin ( angle ) * bn[i];
    angle = ( ( ityp ) ( f2 + f3 ) * theta_ac ) / ( ityp ) ( f1 + f2 + f3 );
    #pragma omp parallel for num_threads(3)
    for ( i = 0; i < 3; ++i )
        ac[i] = cos ( angle ) * a_xyz[i] + sin ( angle ) * cn[i];
    /*
    Determine the angular distance between AB and AC.
    */
    theta_bc = sphere01_distance_xyz ( ab, ac );
    /*
    Polarize AC = ACP + ACN, normalize ACN.
    */
    r8vec_polarize ( 3, ac, ab, acn, acp );
    norm = r8vec_norm ( 3, acn );
    for ( i = 0; i < 3; ++i )
        acn[i] /= norm;
    /*
    The interval between AB and AC is marked by F2+F3+1 vertices 0 through F2+F3.
    */
    angle = ( ( ityp ) ( f3 ) * theta_bc ) / ( ityp ) ( f2 + f3 );

    #pragma omp parallel for num_threads(3)
    for ( i = 0; i < 3; ++i )
        node_xyz[i] = cos ( angle ) * ab[i] + sin ( angle ) * acn[i];
    return node_xyz;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _sphere01_triangle_sample ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE01_TRIANGLE_SAMPLE: sample points from triangle on unit sphere.
  Discussion:
    The sphere has center 0 and radius 1.
    A spherical triangle on the surface of the unit sphere contains those
    points with radius R = 1, bounded by the vertices V1, V2, V3.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    24 September 2010
  Author:
    John Burkardt
  Reference:
    James Arvo,
    Stratified sampling of spherical triangles,
    Computer Graphics Proceedings, Annual Conference Series,
    ACM SIGGRAPH '95, pages 437-438, 1995.
  Parameters:
    Input, int N, the number of points.
    Input, double V1[3], V2[3], V3[3], the XYZ coordinates of
    the vertices of the spherical triangle.
    Input/output, int *SEED, a seed for the random number generator.
    Output, double SPHERE01_TRIANGLE_SAMPLE[3*N], the XYZ coordinates of the
    sample points.
*/
{
	const dt3pitpi * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * v1 = s_data->a1;
	ityp * v2 = s_data->a2;
	ityp * v3 = s_data->a3;
	int * seed = s_data->a4;
	
    ityp a;
    ityp alpha;
    ityp area;
    ityp area_hat;
    ityp b;
    ityp beta;
    ityp c;
    ityp gamma;
    dim_typ i, j;
    ityp norm;
    ityp q;
    ityp s;
    ityp t;
    ityp u;
    ityp v;
    ityp v31[3];
    ityp v4[3];
    ityp v42[3];
    ityp w;
    ityp *x;
    ityp xsi1;
    ityp xsi2;
    ityp z;

    sphere01_triangle_vertices_to_sides ( v1, v2, v3, &a, &b, &c );
    sphere01_triangle_sides_to_angles ( a, b, c, &alpha, &beta, &gamma );
    area = sphere01_triangle_angles_to_area ( alpha, beta, gamma );

    x = ( ityp * ) malloc ( 3 * n * sizeof ( ityp ) );

    for ( j = 0; j < n; ++j)
    {
        /*
        Select the new area.
        */
        xsi1 = r8_uniform_01 ( seed );
        area_hat = xsi1 * area;
        /*
        Compute the sine and cosine of the angle phi.
        */
        s = sin ( area_hat - alpha );
        t = cos ( area_hat - alpha );
        /*
        Compute the pair that determines beta_hat.
        */
        u = t - cos ( alpha );
        v = s + sin ( alpha ) * cos ( c );
        /*
        Q is the cosine of the new edge length b_hat.
        */
        q = ( ( v * t - u * s ) * cos ( alpha ) - v )/ ( ( v * s + u * t ) * sin ( alpha ) );
        /*
        Occasionally we get a Q value out of bounds.
        */
        q = MAX ( q, - 1.00 );
        q = MIN ( q, + 1.00 );
        /*
        V31 = normalized ( V3 - ( V3 dot V1 ) * V1 )
        */
        w = r8vec_dot_product ( 3, v3, v1 );
        #pragma omp parallel for num_threads(3)
        for ( i = 0; i < 3; ++i )
        v31[i] -= w * v1[i];
        norm = r8vec_norm ( 3, v31 );
        #pragma omp parallel for num_threads(3)
        for ( i = 0; i < 3; ++i)
            v31[i] /= norm;
        /*
        V4 is the third vertex of the subtriangle V1, V2, V4.
        */
        #pragma omp parallel for num_threads(3)
        for ( i = 0; i < 3; ++i )
            v4[i] = q * v1[i] + sqrt ( 1.00 - q * q ) * v31[i];
        /*
        Select cos theta, which will sample along the edge from V2 to V4.
        */
        xsi2 = r8_uniform_01 ( seed );
        z = 1.00 - xsi2 * ( 1.00 - r8vec_dot_product ( 3, v4, v2 ) );
        /*
        V42 = normalized ( V4 - ( V4 dot V2 ) * V2 )
        */
        w = r8vec_dot_product ( 3, v4, v2 );
        #pragma omp parallel for num_threads(3)
        for ( i = 0; i < 3; ++i )
            v42[i] = v4[i] - w * v2[i];
        norm = r8vec_norm ( 3, v42 );
        #pragma omp parallel for num_threads(3)
        for ( i = 0; i < 3; ++i )
            v42[i] /= norm;
        /*
        Construct the point.
        */
        #pragma omp parallel for num_threads(3)
        for ( i = 0; i < 3; ++i )
            x[i+j*3] = z * v2[i] + sqrt ( 1.00 - z * z ) * v42[i];
    }
    return x;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _sphere01_triangle_sides_to_angles ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE01_TRIANGLE_SIDES_TO_ANGLES: spherical triangle angles on the unit sphere.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 September 2010
  Author:
    John Burkardt
  Parameters:
    Input, double AS, BS, CS, the (geodesic) length of the sides of the
    triangle.
    Output, double *A, *B, *C, the spherical angles of the triangle.
    Angle A is opposite the side of length AS, and so on.
*/
{
	const _3it3pit * const s_data = data;
	const register ityp as = s_data->a0;
	const register ityp bs = s_data->a1;
	const register ityp cs = s_data->a2;
	ityp * a = s_data->a3;
	ityp * b = s_data->a4;
	ityp * c = s_data->a5;
	
    const register ityp ssu = ( as + bs + cs ) / 2.00;
    *a = 2.00 * atan ( ( sin ( ssu - bs ) * sin ( ssu - cs ) ) /( sin ( ssu ) * sin ( ssu - as )     ) );
    *b = 2.00 * atan ( ( sin ( ssu - as ) * sin ( ssu - cs ) ) /( sin ( ssu ) * sin ( ssu - bs )     ) );
    *c = 2.00 * atan ( sqrt ( ( sin ( ssu - as ) * sin ( ssu - bs ) ) /( sin ( ssu ) * sin ( ssu - cs )     ) ) );
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _sphere01_triangle_vertices_to_angles ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE01_TRIANGLE_VERTICES_TO_ANGLES: angles of a spherical triangle on the unit sphere.
  Discussion:
    A unit sphere centered at 0 in 3D satisfies the equation:
      X*X + Y*Y + Z*Z = 1
    A spherical triangle is specified by three points on the surface
    of the sphere.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 September 2010
  Author:
    John Burkardt
  Parameters:
    Input, double V1[3], V2[3], V3[3], the vertices of the triangle.
    Output, double *A, *B, *C, the angles of the spherical triangle.
*/
{
	ityp ** const a_data = data;
	ityp * v1 = a_data[0];
	ityp * v2 = a_data[1];
	ityp * v3 = a_data[2];
	ityp * a = a_data[3];
	ityp * b = a_data[4];
	ityp * c = a_data[5];
	
    ityp as;
    ityp bs;
    ityp cs;
    /*
    Compute the lengths of the sides of the spherical triangle.
    */
    sphere01_triangle_vertices_to_sides ( v1, v2, v3, &as, &bs, &cs );
    /*
    Get the spherical angles.
    */
    sphere01_triangle_sides_to_angles ( as, bs, cs, a, b, c );

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _sphere01_triangle_vertices_to_area ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE01_TRIANGLE_VERTICES_TO_AREA: area of a spherical triangle on the unit sphere.
  Discussion:
    A sphere centered at 0 in 3D satisfies the equation:
      X*X + Y*Y + Z*Z = 1
    A spherical triangle is specified by three points on the surface
    of the sphere.
    The area formula is known as Girard's formula.
    The area of a spherical triangle on the unit sphere is:
      AREA = A + B + C - M_PI
    where A, B and C are the (surface) angles of the triangle.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 September 2010
  Author:
    John Burkardt
  Parameters:
    Input, double V1[3], V2[3], V3[3], the vertices of the triangle.
    Output, double STRI_VERTICES_TO_AREA_3D, the area of the
    spherical triangle.
*/
{
	static ityp result = MAX_VAL;
	
	ityp ** const a_data = data;
	ityp * v1 = a_data[0];
	ityp * v2 = a_data[1];
	ityp * v3 = a_data[2];
	
    ityp area;
    ityp a;
    ityp as;
    ityp b;
    ityp bs;
    ityp c;
    ityp cs;
    /*
    Compute the lengths of the sides of the spherical triangle.
    */
    sphere01_triangle_vertices_to_sides ( v1, v2, v3, &as, &bs, &cs );
    /*
    Get the spherical angles.
    */
    sphere01_triangle_sides_to_angles ( as, bs, cs, &a, &b, &c );
    /*
    Get the area
    */
    
    result = sphere01_triangle_angles_to_area ( a, b, c );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _sphere01_triangle_vertices_to_centroid ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE01_TRIANGLE_VERTICES_TO_CENTROID: spherical triangle centroid on the unit sphere.
  Discussion:
    A unit sphere centered at 0 in 3D satisfies the equation:
      X*X + Y*Y + Z*Z = 1
    A spherical triangle is specified by three points on the sphere.
    The (true) centroid of a spherical triangle is the point
      VT = (XT,YT,ZT) = Integral ( X, Y, Z ) dArea / Integral 1 dArea
    Note that the true centroid does NOT, in general, lie on the sphere.
    The "flat" centroid VF is the centroid of the planar triangle defined by
    the vertices of the spherical triangle.
    The "spherical" centroid VS of a spherical triangle is computed by
    the intersection of the geodesic bisectors of the triangle angles.
    The spherical centroid lies on the sphere.
    VF, VT and VS lie on a line through the center of the sphere.  We can
    easily calculate VF by averaging the vertices, and from this determine
    VS by normalizing.
 (Of course, we still will not have actually computed VT, which lies
    somewhere between VF and VS!)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 September 2010
  Author:
    John Burkardt
  Parameters:
    Input, double V1[3], V2[3], V3[3], the vertices of the triangle.
    Output, double SPHERE01_TRIANGLE_VERTICES_TO_CENTROID[3], the coordinates of the
    "spherical centroid" of the spherical triangle.
*/
{
	ityp ** const a_data = data;
	ityp * v1 = a_data[0];
	ityp * v2 = a_data[1];
	ityp * v3 = a_data[2];
	
    # define DIM_NUM 3

    dim_typ i;
    ityp norm;
    ityp *vs;

    vs = ( ityp * ) malloc ( 3 * sizeof ( ityp ) );

    #pragma omp parallel for num_threads(3)
    for ( i = 0; i < DIM_NUM; ++i )
        vs[i] = ( v1[i] + v2[i] + v3[i] ) / 3.00;

    norm = r8vec_norm ( DIM_NUM, vs );

    #pragma omp parallel for num_threads(3)
    for ( i = 0; i < DIM_NUM; ++i )
        vs[i] /= norm;

    return vs;
    # undef DIM_NUM
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _sphere01_triangle_vertices_to_midpoints ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE01_TRIANGLE_VERTICES_TO_MIDPOINTS gets the midsides of a spherical triangle.
  Discussion:
    The points are assumed to lie on the unit sphere.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    23 September 2010
  Author:
    John Burkardt
  Parameters:
    Input, double V1[3], V2[3], V3[3], the vertices of the triangle.
    Output, double M1[3], M2[3], M3[3], the coordinates of
    the midpoints of the sides of the spherical triangle.
*/
{
	ityp ** const a_data = data;
	ityp * v1 = a_data[0];
	ityp * v2 = a_data[1];
	ityp * v3 = a_data[2];
	ityp * m1 = a_data[3];
	ityp * m2 = a_data[4];
	ityp * m3 = a_data[5];
	
    dim_typ i;
    ityp norm;

    #pragma omp parallel for num_threads(3)
    for ( i = 0; i < 3; ++i )
        m1[i] = ( v1[i] + v2[i] ) / 2.00;
    norm = r8vec_norm ( 3, m1 );
    #pragma omp parallel for num_threads(3)
    for ( i = 0; i < 3; ++i )
    {
        m1[i] /= norm;
        m2[i] = ( v2[i] + v3[i] ) / 2.00;
    }
    norm = r8vec_norm ( 3, m2 );
    #pragma omp parallel for num_threads(3)
    for ( i = 0; i < 3; ++i )
    {
        m2[i] /= norm;
        m3[i] = ( v3[i] + v1[i] ) / 2.00;
    }
    norm = r8vec_norm ( 3, m3 );
    for ( i = 0; i < 3; ++i )
        m3[i] /= norm;

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _sphere01_triangle_vertices_to_sides ( void * data)
/******************************************************************************/
/*
  Purpose:
    SPHERE01_TRIANGLE_VERTICES_TO_SIDES_3D: spherical triangle sides on the unit sphere.
  Discussion:
    We can use the ACOS system call here, but the acos routine
    will automatically take care of cases where the input argument is
 (usually slightly) out of bounds.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    25 September 2010
  Author:
    John Burkardt
  Parameters:
    Input, double V1[3], V2[3], V3[3], the vertices of the spherical
    triangle.
    Output, double *AS, *BS, *CS, the (geodesic) length of the sides of the
    triangle.
*/
{
	ityp ** const a_data = data;
	ityp * v1 = a_data[0];
	ityp * v2 = a_data[1];
	ityp * v3 = a_data[2];
	ityp * as = a_data[3];
	ityp * bs = a_data[4];
	ityp * cs = a_data[5];
	
    *as = acos( r8vec_dot_product ( 3, v2, v3 ) );
    *bs = acos ( r8vec_dot_product ( 3, v3, v1 ) );
    *cs = acos ( r8vec_dot_product ( 3, v1, v2 ) );
    return NULL;
}

#endif
