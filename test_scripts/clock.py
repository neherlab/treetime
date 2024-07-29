from Bio import AlignIO
import numpy as np
from lib import Graph, graph_from_nwk_str, GraphNodeBackward, GraphNodeForward, Node, ClockSet, ClockModel
from typing import List, Dict
from payload import NodePayload, EdgePayload

def propagate_averages(Q_src: ClockSet, branch_value: float, branch_variance: float) -> ClockSet:
  denom = 1.0/(1+branch_variance*Q_src.norm)
  # denom = 1/q/(1/q + bv)
  # q = sum_c 1/q_c/(1/q_c + bv_c)
  # 1/q = 1/(sum_c 1/q_c/(1/q_c + bv_c)) =
  return ClockSet(
    t_sum  = Q_src.t_sum*denom,  # Eq 11 in Neher 2018 -- contribution of children doesn't change
    tsq_sum   = Q_src.tsq_sum - branch_variance*Q_src.t_sum**2*denom, # Eq. 13
    d_sum  = (Q_src.d_sum + branch_value*Q_src.norm)*denom, # Eq. 12 -- add branch_value norm times
    dt_sum = (Q_src.dt_sum + branch_value*Q_src.t_sum - branch_variance*Q_src.t_sum*(Q_src.d_sum+Q_src.norm*branch_value))*denom, # Eq. 14
    # Eq. A.2
    dsq_sum = (Q_src.dsq_sum + 2*branch_value*Q_src.d_sum + branch_value**2*Q_src.norm \
               - branch_variance*(Q_src.d_sum**2 + 2*branch_value*Q_src.d_sum*Q_src.norm + branch_value**2*Q_src.norm**2)*denom),
    norm   = Q_src.norm*denom
  )


def clock_regression_backward(graph: Graph) -> None:
  bv = lambda e: 0.0
  def clock_backward(n: GraphNodeBackward) -> None:
    n.payload.clock.from_children = {}
    if n.is_leaf:
      t = n.payload.clock.date
      n.payload.clock.to_parent = ClockSet(t_sum=t, tsq_sum=t**2, d_sum=0, dsq_sum=0, dt_sum=0, norm=1.0)
    else:
      Q_dest = ClockSet()
      for c,e in n.children:
        n.payload.clock.from_children[c.name] = propagate_averages(c.clock.to_parent, branch_value=e.branch_length,
                                                                            branch_variance=bv(e))
        Q_dest.add(n.payload.clock.from_children[c.name])

      n.payload.clock.to_parent=Q_dest

  graph.par_iter_backward(clock_backward)

def clock_regression_forward(graph: Graph) -> None:
  bv = lambda e: 0.0
  def clock_forward(n: GraphNodeForward) -> None:
    Q = n.payload.clock.to_parent
    Q_tot = ClockSet(t_sum=Q.t_sum, tsq_sum=Q.tsq_sum, d_sum=Q.d_sum, dsq_sum=Q.dsq_sum, dt_sum=Q.dt_sum, norm=Q.norm)
    if not n.is_root:
      p,e = n.parents[0]
      Q_tot.add(propagate_averages(p.clock.to_children[n.payload.name], branch_value=e.branch_length, branch_variance=bv(e)))

    n.payload.clock.to_children = {}
    print(n.payload.clock.from_children)
    for c, Q in n.payload.clock.from_children.items():
      Q_to_children = ClockSet(t_sum=Q_tot.t_sum, tsq_sum=Q_tot.tsq_sum, d_sum=Q_tot.d_sum,
                               dsq_sum=Q_tot.dsq_sum, dt_sum=Q_tot.dt_sum, norm=Q_tot.norm)

      Q_to_children.subtract(Q)
      n.payload.clock.to_children[c] = Q_to_children

    n.payload.clock.total = Q_tot

  graph.par_iter_forward(clock_forward)


def clock_model(Q: ClockSet) -> ClockModel:
  det = Q.tsq_sum*Q.norm - Q.t_sum**2
  if det>0: # make sure there is variation in dates
    rate = (Q.dt_sum*Q.norm - Q.t_sum*Q.d_sum) / det
  else:
    raise ValueError("No variation in sampling dates! Please specify your clock rate explicitly.")

  intercept = (Q.t_sum - Q.t_sum*rate)/Q.norm
  estimator_hessian = np.array([[Q.tsq_sum, Q.t_sum], [Q.t_sum, Q.norm]])
  chisq = 0.5*(Q.dsq_sum*Q.norm - Q.d_sum**2 - (Q.dt_sum*Q.norm - Q.d_sum*Q.t_sum)**2/(Q.tsq_sum*Q.norm - Q.t_sum**2))/Q.norm
  return ClockModel(rate=rate, intercept=intercept, hessian=estimator_hessian, chisq=chisq)

def clock_model_fixed_rate(Q: ClockSet, rate: float) -> ClockModel:

  intercept = (Q.t_sum - Q.t_sum*rate)/Q.norm
  estimator_hessian = np.array([[Q.tsq_sum, Q.t_sum], [Q.t_sum, Q.norm]])

  return ClockModel(rate=rate, intercept=intercept, hessian=estimator_hessian)



if __name__=="__main__":
  fname_nwk = 'test_scripts/data/tree.nwk'
  fname_seq = 'test_scripts/data/sequences.fasta'
#   fname_metadata = 'test_scripts/data/metadata.tsv'
  # ((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;

  dates = {"A":2013, "B": 2022, "C":2017, "D": 2005}
  div = {"A":0.2, "B": 0.3, "C":0.25, "D": 0.17}
  t = np.sum([x for x in dates.values()])
  tsq = np.sum([x**2 for x in dates.values()])
  dt = np.sum([div[c]*dates[c] for c in div])
  d = np.sum([div[c] for c in div])
  print((dt*4-d*t)/(tsq*4-t**2))

  with open(fname_nwk) as fh:
    nwkstr = fh.read()
  G = graph_from_nwk_str(nwk_string=nwkstr, node_payload_factory=NodePayload, edge_payload_factory=EdgePayload)

  for n in G.get_leaves():
    n.payload().clock.date = dates[n.payload().name]

  clock_regression_backward(G)
  clock_regression_forward(G)

  #
  root = G.get_one_root()
  print(root.payload().clock)

  print(clock_model(root.payload().clock.total))