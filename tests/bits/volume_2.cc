//----------------------------  volume_2.cc  ---------------------------
//    volume_2.cc,v 1.1 2003/10/19 22:29:39 wolf Exp
//    Version: 
//
//    Copyright (C) 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  volume_2.cc  ---------------------------


// like a union of the normals_* and q_point_sum_* tests: integrating
// \vec x times the normal over the surface of any body should yield
// the volume of this body times the space dimension
//
// in contrast to volume_1, do this with a Q1, Q2, and Q4 mapping
// associated with the geometry
//
// this used to fail in 3d, since there were assumptions on face
// orientations present that had to be weeded out

#include "../tests.h"
#include <base/logstream.h>
#include <base/quadrature_lib.h>
#include <grid/tria.h>
#include <grid/tria_boundary_lib.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_generator.h>
#include <dofs/dof_handler.h>
#include <fe/fe_q.h>
#include <fe/fe_values.h>
#include <fe/mapping_q.h>

#include <fstream>



template <int dim>
void check (const Triangulation<dim> &tria,
            const unsigned int        order)
{
  MappingQ<dim> mapping(order);
  FE_Q<dim> fe(1);
  DoFHandler<dim> dof_handler (tria);
  dof_handler.distribute_dofs (fe);

  QGauss3<dim-1> q_face;
  
  FEFaceValues<dim>    fe_face_values (mapping, fe, q_face,
                                       update_normal_vectors |
                                       update_q_points |
                                       update_JxW_values);
  FESubfaceValues<dim> fe_subface_values (mapping, fe, q_face,
                                          update_normal_vectors |
                                          update_q_points |
                                          update_JxW_values);

  double v1=0, v2=0;
  for (typename DoFHandler<dim>::active_cell_iterator
         cell = dof_handler.begin_active();
       cell!=dof_handler.end(); ++cell)
    {

                                       // first integrate over faces
                                       // and make sure that the
                                       // result of the integration is
                                       // close to zero
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
        if (cell->at_boundary(f))
          {
            fe_face_values.reinit (cell, f);
            for (unsigned int q=0; q<q_face.n_quadrature_points; ++q)
              v1 += (fe_face_values.normal_vector(q) *
                     fe_face_values.quadrature_point(q)) *
                    fe_face_values.JxW(q);
          }

                                       // now same for subface
                                       // integration
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
        if (cell->at_boundary(f))
          for (unsigned int sf=0; sf<GeometryInfo<dim>::subfaces_per_face; ++sf)
            {
              fe_subface_values.reinit (cell, f, sf);
              for (unsigned int q=0; q<q_face.n_quadrature_points; ++q)
                v2 += (fe_subface_values.normal_vector(q) *
                       fe_subface_values.quadrature_point(q)) *
                      fe_subface_values.JxW(q);
            }
    }

  Assert (std::fabs(v1-v2)/v1 < 4e-4, ExcInternalError());
  deallog << " face integration: "
          << v1 / dim
          << std::endl;
  deallog << " subface integration: "
          << v2 / dim
          << std::endl;
}


int main () 
{
  std::ofstream logfile("volume_2.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  {  
    Triangulation<2> coarse_grid;
    GridGenerator::hyper_ball (coarse_grid);
    static const HyperBallBoundary<2> boundary;
    coarse_grid.set_boundary (0, boundary);
    check (coarse_grid, 1);
    check (coarse_grid, 2);
    check (coarse_grid, 4);
  }
  
  {  
    Triangulation<3> coarse_grid;
    GridGenerator::hyper_ball (coarse_grid);
    static const HyperBallBoundary<3> boundary;
    coarse_grid.set_boundary (0, boundary);
    check (coarse_grid, 1);
    check (coarse_grid, 2);
    check (coarse_grid, 3);
  }
  
}

  
  
