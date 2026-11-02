"""Display-only ARF grouping and probability-preserving smoothing."""
import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'diagnostics/arf_display'))
from combined import category
from main_figure import smooth_simplex
class DisplayContracts(unittest.TestCase):
 def test_all_source_states(self):
  expected={('H0','C0'):3,('H0','C1'):1,('H0','CU'):3,('H1','C0'):0,('H1','C1'):2,('H1','CU'):0,('HU','C0'):3,('HU','C1'):1,('HU','CU'):3}
  self.assertTrue(all(category(*k)==v for k,v in expected.items()))
 def test_smoothing_simplex(self):
  raw=np.random.default_rng(72).dirichlet(np.ones(4),size=72).T;before=raw.copy();x,z=smooth_simplex(raw)
  self.assertTrue(np.array_equal(raw,before));self.assertTrue(np.all(z>=0));self.assertTrue(np.allclose(z.sum(axis=0),1));self.assertEqual((x[0],x[-1]),(0,72))
 def test_constant_unknown_and_zero_states(self):
  raw=np.zeros((4,72));raw[3]=1;_,z=smooth_simplex(raw)
  self.assertTrue(np.all(z[3]==1));self.assertTrue(np.all(z[:3]==0))
 def test_edge_normalization(self):
  raw=np.zeros((4,72));raw[0,0]=1;raw[3,1:]=1;_,z=smooth_simplex(raw)
  self.assertAlmostEqual(z[0,0],1/3);self.assertTrue(np.allclose(z.sum(axis=0),1))
if __name__=='__main__':unittest.main()
