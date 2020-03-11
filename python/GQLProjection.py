import numpy as np

from dedalus.core.field import Operand
from dedalus.core.operators import Operator, FutureField, Interpolate


class GQLProjection(Operator, FutureField):
    """
    Projection operator for generalized quasilinear approximation
    
    
    """
    def __init__(self, arg, cutoff, subspace, dim=None,**kw):
        arg = Operand.cast(arg)
        super().__init__(arg,**kw)
        self.cutoff = cutoff

        
        # by default, execute GQL on all but the last dimension
        if not dim:
            self.dim = self.domain.dim - 1
        else:
            self.dim = dim
        
        local_coeff = self.domain.all_elements()
        low_mask = np.ones(self.domain.local_coeff_shape, dtype='bool')

        for i in range(self.dim):
            low_mask &= (np.abs(local_coeff[i]) <= cutoff[i])
        if subspace == 'high' or subspace == 'h':
            self.mask = ~low_mask
            subspace_name = 'h'
        elif subspace == 'low' or subspace == 'l':
            self.mask = low_mask
            subspace_name = 'l'
        else:
            raise ValueError("Subspace must be high/h or low/l, not {}".format(subspace))

        cutoff_str = ",".join([str(i) for i in cutoff]) 
        self.name = 'Proj[({}),{}]'.format(cutoff_str,subspace_name)

    def meta_constant(self, axis):
        # Preserve constancy
        return self.args[0].meta[axis]['constant']

    def check_conditions(self):
        """Projection must be in coefficient space""" 
        return self.args[0].layout is self._coeff_layout
    
    def operate(self, out):
        for i in range(self.dim):
            self.args[0].require_layout('c')
        out.data[:] = self.args[0].data
        out.data *= self.mask
