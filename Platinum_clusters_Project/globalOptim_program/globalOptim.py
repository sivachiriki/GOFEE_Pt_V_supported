import numpy as np
from time import time

from ase.io import read, write, Trajectory
from ase.calculators.singlepoint import SinglePointCalculator

from globalOptim_baseClass import relaxGPAW, singleGPAW, relax_VarianceBreak, globalOptim_baseClass

import os

import ase.parallel as mpi
world = mpi.world


class globalOptim(globalOptim_baseClass):

    def __init__(self, traj_namebase, MLmodel, startGenerator, mutationSelector, calc=None, startStructures=None, slab=None, population_size=5, kappa=2, Niter=50, Ninit=2, Ncand_min=30, max_similarity_dist=0.05, min_certainty=0.8, zlim=None, dualPoint=False, errorRelax=False, relaxFinalPop=False, Nuse_pop_as_candidates=1):

        globalOptim_baseClass.__init__(self,
                                       traj_namebase=traj_namebase,
                                       MLmodel=MLmodel,
                                       startGenerator=startGenerator,
                                       mutationSelector=mutationSelector,
                                       calc=calc,
                                       startStructures=startStructures,
                                       slab=slab,
                                       population_size=population_size,
                                       kappa=kappa,
                                       Niter=Niter,
                                       Ninit=Ninit,
                                       Ncand_min=Ncand_min,
                                       max_similarity_dist=max_similarity_dist,
                                       min_certainty=min_certainty,
                                       zlim=zlim,
                                       dualPoint=dualPoint,
                                       errorRelax=errorRelax,
                                       relaxFinalPop=relaxFinalPop,
                                       Nuse_pop_as_candidates=Nuse_pop_as_candidates)

    def relaxTrue(self, a):
        #h=5+'s'
        #print('i came till top',flush=True)
        pos = a.positions
        if self.master:
            pos = a.positions
         #   print('i came till bet',flush=True)
        self.comm.broadcast(pos, 0)
        a.positions = pos
        #print('i came till belo',flush=True)
        self.comm.barrier()

        # Relax
        label = self.traj_namebase + '{}'.format(self.traj_counter)
        #print('i came till',flush=True)
        a_relaxed = relaxGPAW(a, label, calc=self.calc)
        #print('relaxed',flush=True)
        # Add sampled trajectory to training data.
        self.add_trajectory_to_training(label+'_lcao.traj')

        self.traj_counter += 1
        Erelaxed = a_relaxed.get_potential_energy()
        Frelaxed = a_relaxed.get_forces()
        return a_relaxed, Erelaxed, Frelaxed
            
    def singlePoint(self, anew):
        a = anew.copy()
        
        # Save structure with ML-energy
        if self.master:
            self.writer_spPredict.write(a)

        # broadcast structure, so all cores have the same
        pos = a.positions
        if self.master:
            pos = a.positions
        self.comm.broadcast(pos, 0)
        a.positions = pos
        self.comm.barrier()

        # Perform single-point
        label =  self.traj_namebase + '{}'.format(self.traj_counter)
        E, F = singleGPAW(a, label, calc=self.calc)
        self.comm.barrier()

        # save structure for training
        results = {'energy': E, 'forces': F}
        calc_sp = SinglePointCalculator(a, **results)
        a.set_calculator(calc_sp)
        self.a_add.append(a)

        # Save to spTrain-trajectory
        self.writer_spTrain.write(a, energy=E, forces=F)

        self.traj_counter += 1
        return E, F

    
    
