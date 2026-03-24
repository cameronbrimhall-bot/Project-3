#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 14:34:19 2021

@author: kendrick
"""

import numpy as np

# compute unknown displacements 
def ComputeDisplacements(K, F, n_unknowns):
    # extract submatrix of unknowns
    K11 = K[0:n_unknowns,0:n_unknowns]
    F1 = F[0:n_unknowns]
    
    d = np.linalg.solve(K11,F1)
    
    return d

# postprocess the forces at known displacement nodes
def PostprocessReactions(K, d, F, n_unknowns, nodes):
    # These are computed net forces and do not
    # take into account external loads applied
    # at these nodes
    F = np.matmul(K[n_unknowns:,0:n_unknowns], d)
    
    # Postprocess the reactions
    for node in nodes:
        if node.xidx >= n_unknowns:
            node.AddReactionXForce(F[node.xidx-n_unknowns][0] - node.xforce_external)
        if node.yidx >= n_unknowns:
            node.AddReactionYForce(F[node.yidx-n_unknowns][0] - node.yforce_external)
        
    return F

# determine internal member loads
def ComputeMemberForces(bars):
    """
    Computes internal axial member loads for each bar using Equation 14-23.
    Iterates through each member, extracts properties and node displacements,
    and performs the matrix calculation to solve for the axial load.
    """
    for bar in bars:
        # ii. Extract properties
        E = bar.E
        A = bar.A
        L = bar.Length()
        [lx, ly] = bar.LambdaTerms()
        
        # iii. Extract displacements from near (init) and far (end) nodes
        # d = [u_near, v_near, u_far, v_far]
        d = np.array([
            bar.init_node.xdisp,
            bar.init_node.ydisp,
            bar.end_node.xdisp,
            bar.end_node.ydisp
        ])
        
        # iv. Perform calculation based on Equation 14-23:
        # q = (AE/L) * [-lx, -ly, lx, ly] * {d}
        transformation_vector = np.array([-lx, -ly, lx, ly])
        axial_force = (A * E / L) * np.dot(transformation_vector, d)
        
        # v. Store the axial load
        bar.axial_load = axial_force
        pass
    
# compute the normal stresses
def ComputeNormalStresses(bars):
   for bar in bars:
        if bar.A != 0:
            bar.normal_stress = bar.axial_load / bar.A
        else:
            bar.normal_stress = 0.0
        pass

# compute the critical buckling load of a member
def ComputeBucklingLoad(bars):
  """
    Computes the critical buckling load P_cr = (pi^2 * E * I) / (KL)^2.
    For trusses, K = 1.0.
    """
  for bar in bars:
        E = bar.E
        # Use the weak axis (Iu) for the most conservative buckling calculation
        I = bar.Iu 
        L = bar.Length()
        K = 1.0
        
        # P_cr calculation
        critical_load = (np.pi**2 * E * I) / (K * L)**2
        bar.buckling_load = critical_load
        pass
