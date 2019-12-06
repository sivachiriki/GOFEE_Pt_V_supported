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
        pos = a.positions
        if self.master:
            pos = a.positions
        self.comm.broadcast(pos, 0)
        a.positions = pos
        self.comm.barrier()

        label = self.traj_namebase + '{}'.format(self.traj_counter)
        
        pos_relaxed = a.positions
        Erelaxed = np.zeros(1)
        Frelaxed = np.zeros((self.Natoms,3))
        if self.master:
            # Relax
            a_relaxed = relaxGPAW(a, label, calc=self.calc)
            pos_relaxed = a_relaxed.positions
            Erelaxed = np.array([a_relaxed.get_potential_energy()])
            Frelaxed = a_relaxed.get_forces()
        self.comm.broadcast(pos_relaxed, 0)
        self.comm.broadcast(Erelaxed, 0)
        Erelaxed = Erelaxed[0]
        self.comm.broadcast(Frelaxed, 0)
        a_relaxed = a.copy()
        a_relaxed.positions = pos_relaxed        
        self.comm.barrier()

        # Add sampled trajectory to training data.
        self.add_trajectory_to_training(label+'_lcao.traj')

        self.traj_counter += 1
        return a_relaxed, Erelaxed, Frelaxed
    
    def singlePoint(self, anew):
        a = anew.copy()
        
        # Save structure with ML-energy
        if self.master:
            self.writer_spPredict.write(a)
        #self.writer_spPredict.write(a)
        
        # broadcast structure, so all cores have the same
        pos = a.positions
        if self.master:
            pos = a.positions
        self.comm.broadcast(pos, 0)
        a.positions = pos
        self.comm.barrier()

        
        # Perform single-point - With error handling
        E = np.zeros(1)
        F = np.zeros((self.Natoms,3))
        try:
            error_int = np.array([0]).astype(int)
            if self.master:
                label =  self.traj_namebase + '{}'.format(self.traj_counter)
                E, F = singleGPAW(a, label, calc=self.calc)
                E = np.array([E])
        except Exception as runtime_err:
            error_int = np.array([1]).astype(int)
            print('Error cought:', runtime_err)

        world.barrier()
        self.comm.broadcast(error_int, 0)
        if error_int[0]==1:
            raise RuntimeError('DFTB failed')
        world.barrier()
        
        self.comm.broadcast(E, 0)
        E = E[0]
        self.comm.broadcast(F, 0)
        
        
        # save structure for training
        a.energy = E
        results = {'energy': E}
        calc = SinglePointCalculator(a, **results)
        a.set_calculator(calc)
        self.a_add.append(a)

        # Save to spTrain-trajectory
        self.writer_spTrain.write(a, energy=E, forces=F)

        self.traj_counter += 1
        return E, F
