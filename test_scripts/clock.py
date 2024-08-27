import numpy as np
from lib import Graph, graph_from_nwk_str, GraphNodeBackward, GraphNodeForward, Node, ClockSet, ClockModel, FindRootResult, GraphEdgeKey
from typing import List, Dict
from payload import NodePayload, EdgePayload

def propagate_averages(Q_src: ClockSet, branch_value: float, branch_variance: float) -> ClockSet:
  denom = 1.0/(1+branch_variance*Q_src.norm)
  return ClockSet(
    t_sum  = Q_src.t_sum*denom,  # Eq 11 in Neher 2018 -- contribution of children doesn't change
    tsq_sum   = Q_src.tsq_sum - branch_variance*Q_src.t_sum**2*denom, # Eq. 13
    d_sum  = (Q_src.d_sum + branch_value*Q_src.norm)*denom, # Eq. 12 -- add branch_value norm times
    dt_sum = Q_src.dt_sum + branch_value*Q_src.t_sum - branch_variance*Q_src.t_sum*(Q_src.d_sum+Q_src.norm*branch_value)*denom, # Eq. 14
    # Eq. A.2
    dsq_sum = (Q_src.dsq_sum + 2*branch_value*Q_src.d_sum + branch_value**2*Q_src.norm \
               - branch_variance*(Q_src.d_sum**2 + 2*branch_value*Q_src.d_sum*Q_src.norm + branch_value**2*Q_src.norm**2)*denom),
    norm   = Q_src.norm*denom
  )

def leaf_contribution(tip_value: float) -> ClockSet:
  return ClockSet(
    t_sum  = tip_value,
    tsq_sum = tip_value**2,
    d_sum  = 0,
    dt_sum = 0,
    dsq_sum = 0,
    norm = 1
  )


def clock_regression_backward(graph: Graph) -> None:
  bv = lambda e: graph.clock_mode['variance_factor']*e.branch_length +  graph.clock_mode['variance_offset'] # should be configurable
  def clock_backward(n: GraphNodeBackward) -> None:
    n.payload.clock.from_children = {}
    if n.is_leaf:
      Q_dest = leaf_contribution(tip_value=n.payload.clock.date)
    else:
      Q_dest = ClockSet()
      for c,e in n.children:
        Q_from_child = propagate_averages(c.clock.to_parent, branch_value=e.branch_length,
                                                       branch_variance=bv(e))
        Q_dest.add(Q_from_child)
        n.payload.clock.from_children[c.name] = Q_from_child

    n.payload.clock.to_parent=Q_dest
  graph.par_iter_backward(clock_backward)

def clock_regression_forward(graph: Graph) -> None:
  bv = lambda e: graph.clock_mode['variance_factor']*e.branch_length +  graph.clock_mode['variance_offset']
  def clock_forward(n: GraphNodeForward) -> None:
    Q = n.payload.clock.to_parent
    Q_tot = ClockSet(t_sum=Q.t_sum, tsq_sum=Q.tsq_sum, d_sum=Q.d_sum, dsq_sum=Q.dsq_sum, dt_sum=Q.dt_sum, norm=Q.norm)
    if not n.is_root:
      p,e = n.parents[0]
      Q_tot.add(propagate_averages(p.clock.to_children[n.payload.name], branch_value=e.branch_length, branch_variance=bv(e)))

    n.payload.clock.to_children = {}
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

  intercept = (Q.d_sum - Q.t_sum*rate)/Q.norm
  estimator_hessian = np.array([[Q.tsq_sum, Q.t_sum], [Q.t_sum, Q.norm]])
  chisq = 0.5*(Q.dsq_sum*Q.norm - Q.d_sum**2 - (Q.dt_sum*Q.norm - Q.d_sum*Q.t_sum)**2/(Q.tsq_sum*Q.norm - Q.t_sum**2))/Q.norm
  return ClockModel(rate=rate, intercept=intercept, hessian=estimator_hessian, chisq=chisq)

def clock_model_fixed_rate(Q: ClockSet, rate: float) -> ClockModel:

  intercept = (Q.t_sum - Q.t_sum*rate)/Q.norm
  estimator_hessian = np.array([[Q.tsq_sum, Q.t_sum], [Q.t_sum, Q.norm]])

  return ClockModel(rate=rate, intercept=intercept, hessian=estimator_hessian)


def find_best_split(graph: Graph, edge: GraphEdgeKey) -> FindRootResult:
  # get clock data from both ends of the edge
  n = graph.get_node(edge.target())
  n_clock = n.payload().clock.to_parent
  p_clock = graph.get_node(edge.source()).payload().clock.to_children[n.payload().name]

  # precalculate branch values for the edge.
  branch_value = edge.payload().branch_length
  branch_variance =  graph.clock_mode['variance_factor']*branch_value +  graph.clock_mode['variance_offset'] # should be configurable

  # interogate different positions along the branch
  best_chisq = np.inf
  best_x = np.nan
  best_clock = None
  for x in np.linspace(0,1,11): # arbitrary choice for now, should optimize
    Q_tmp = propagate_averages(p_clock, branch_value=branch_value*x, branch_variance=branch_variance*x)
    Q_tmp.add(propagate_averages(n_clock, branch_value=branch_value*(1-x), branch_variance=branch_variance*(1-x)))
    clock = clock_model(Q_tmp)

    if clock.chisq<best_chisq:
      best_chisq=clock.chisq
      best_x=x
      best_clock = clock

  return FindRootResult(edge=edge, split=best_x, clock=best_clock)


def find_best_root(graph: Graph) -> FindRootResult:
  '''
  loop over all nodes, pick the one with the lowest chisq, then optimize position along surrounding brances
  '''
  clock_regression_backward(graph)
  clock_regression_forward(graph)
  root = graph.get_one_root()
  best_root_node = root
  best_chisq = clock_model(root.payload().clock.total).chisq
  best_root = FindRootResult(edge=None, split=0.0, clock=clock_model(root.payload().clock.total))

  # find best node
  for n in graph.get_nodes():
    tmp_chisq = clock_model(n.payload().clock.total).chisq
    if tmp_chisq<best_chisq:
      best_root_node=n
      best_chisq=tmp_chisq

  # check if somehwere on parent branch is better
  if best_root_node!=root:
    res = find_best_split(graph, graph.get_edge(best_root_node.inbound()[0]))
    if res.clock.chisq < best_chisq:
      best_chisq=res.clock.chisq
      best_root = res

  # check if somehwere on a child branch is better
  for e in best_root_node.outbound():
    res = find_best_split(graph, graph.get_edge(e))
    if res.clock.chisq < best_chisq:
      best_chisq=res.clock.chisq
      best_root = res

  return best_root


def reroot_on_best_edge(graph: Graph, best_root) -> None:
  n = graph.get_node(best_root.edge.target())
  p = graph.get_node(best_root.edge.source())

  if best_root.split==0.0:
    new_root = n
  elif best_root.split==1.0:
    new_root = p
  else:
    # create new node that will be the new root
    # set this node as the target of the edge with branch length best_root.split*edge.branch_length
    # set this node as the source of a new edge with branch length (1-best_root.split)*edge.branch_length and target n

  # I am not sure this function exists in the python version.
  # graph.reroot(n, new_root, old_root)


def tests():
  dates = {"A":2013, "B": 2022, "C":2017, "D": 2005}
  div = {"A":0.2, "B": 0.3, "C":0.25, "D": 0.17}

  t = np.sum([x for x in dates.values()])
  tsq = np.sum([x**2 for x in dates.values()])
  dt = np.sum([div[c]*dates[c] for c in div])
  d = np.sum([div[c] for c in div])
  naive_rate = (dt*4-d*t)/(tsq*4-t**2)

  tree = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;"
  G = graph_from_nwk_str(nwk_string=tree, node_payload_factory=NodePayload, edge_payload_factory=EdgePayload)
  for n in G.get_leaves():
    n.payload().clock.date = dates[n.payload().name]

  G.clock_mode = {'variance_factor': 0.0, 'variance_offset':0.0}
  clock_regression_backward(G)
  clock_regression_forward(G)

  root = G.get_one_root()
  clock = clock_model(root.payload().clock.total)

  assert naive_rate==clock.rate

  G.clock_mode = {'variance_factor': 1.0, 'variance_offset':0.0}
  res = find_best_root(G)
  assert res.clock.rate == 0.008095476518345305



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
  print("Naive rate:", (dt*4-d*t)/(tsq*4-t**2))

  with open(fname_nwk) as fh:
    nwkstr = fh.read()
  G = graph_from_nwk_str(nwk_string=nwkstr, node_payload_factory=NodePayload, edge_payload_factory=EdgePayload)
  G.clock_mode = {'variance_factor': 1.0, 'variance_offset':0.0}

  for n in G.get_leaves():
    n.payload().clock.date = dates[n.payload().name]

  clock_regression_backward(G)
  clock_regression_forward(G)

  #
  root = G.get_one_root()
  print(root.payload().clock)

  print(clock_model(root.payload().clock.total))

  res = find_best_root(G)
  print(res)