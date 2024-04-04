# This Python file uses the following encoding: utf-8

# if __name__ == "__main__":
#     pass

from scripts.python.graph import DiGraph
from scripts.python.con_types import FlatConTypes

class Model:
  """
  An optimization model with conversion
  graph
  """

  def __init__(self):
    self._graph = DiGraph()    ## Underlyng graph

    self._vars = []            ## Pointers to various parts of the graph
    self._dvars = []
    self._cons_NL_all = []
    self._cons_NL = {        ## NL + SOS
      "All" : [],
      "Nonlinear" : [],
      "Linear" : [],
      "Logical": [],
      "SOS1": [],
      "SOS2": []}
    self._cons_Flat = {}
    self._cons_Flat_Group = {}
    self._objs_NL = []
    self._objs = []

  def UpdateVar(self, idx, data):
    self._updateNodeData(self._vars, idx, data)

  def UpdateDefVar(self, idx, data):
    self._updateNodeData(self._dvars, idx, data)

  def UpdateNLObj(self, idx, data):
    self._updateNodeData(self._objs_NL, idx, data)

  def UpdateFlatObj(self, idx, data):
    self._updateNodeData(self._objs, idx, data)

  def UpdateNLCon(self, type, idx, data):
    if "nonlin"==type or "lin"==type or "logical"==type:
      assert idx==len(self._cons_NL_all)        ## common list of NL constraints
      data = self._updateNodeData(self._cons_NL_all, idx, data)
      assert "node_index" in data               ## Already added to graph
      if "nonlin"==type:
        self._cons_NL["Nonlinear"].append(data)   ## these just 1x
      elif "lin"==type:
        self._cons_NL["Linear"].append(data)
      else:
        self._cons_NL["Logical"].append(data)
    elif "_sos1"==type:                         ## NL SOS constraints via suffixes
      assert "node_index" in data               ## Already added to graph
      self._cons_NL["SOS1"].append(data)
    elif "_sos2"==type:
      assert "node_index" in data
      self._cons_NL["SOS2"].append(data)
    else:
      raise Exception("Unknown NL constraint type: "+type)

  def UpdateFlatConGroup(self, type, data):
    self._cons_Flat_Group[type] = data

  def UpdateFlatCon(self, type, idx, data):
    if type not in self._cons_Flat:
      self._cons_Flat[type] = []
    data1 = self._updateNodeData(self._cons_Flat[type], idx, data)
    if 0==data1["depth"] \
        and type.startswith('_sos') \
        and "printed" in data1:   ## we need the final status
      self.UpdateNLCon(type, 0, data1)

  def _updateNodeData(self, specnodecnt, idx, data):
    data1, upd = self._updateItemData(specnodecnt, idx, data)
    if (upd):
      idx_g = self._graph.AddNode(data1)
      data1["node_index"] = idx_g
    return data1

  def _updateItemData(self, specnodecnt, idx, data):
    if len(specnodecnt)<=idx:
      specnodecnt.insert(idx, {})
    if (specnodecnt[idx] is None):    ## No such item
      specnodecnt[idx] = {}
    ifEmpty = 0==len(specnodecnt[idx])
    self._updateMap(specnodecnt[idx], data)
    return specnodecnt[idx], ifEmpty

  def _updateMap(self, data1, data2):
    data1.update(data2)

  def AddLinks(self, chunk):
    src_nodes = chunk["src_nodes"]
    dest_nodes = chunk["dest_nodes"]
    for s in src_nodes:
      for d in dest_nodes:
        for sk, svL in s.items():                 ## Actually just 1 key-value pair
          for dk, dvL in d.items():
            if not dk.startswith("dest_cons(") \
                and "src_vars()"!=sk:  ## No links from NL vars. TODO extract group index
              ifsV = int==type(svL)                  ## Scalar
              ifdV = int==type(dvL)
              if ifsV and ifdV:
                self.AddLink(sk, svL, dk, dvL)
              elif ifsV:
                for dv in dvL:
                  self.AddLink(sk, svL, dk, dv)
              elif ifdV:
                for sv in svL:
                  self.AddLink(sk, sv, dk, dvL)
              else:
                assert "CopyLink"==chunk["link_type"] and len(svL)==len(dvL)
                for i in range(len(svL)):
                  self.AddLink(sk, svL[i], dk, dvL[i])

  def AddLink(self, key1, i1, key2, i2):
    si = self.GetLinkNode(key1, i1)
    di = self.GetLinkNode(key2, i2)
    self._graph.AddLink(si["node_index"], di["node_index"])

  def GetLinkNode(self, type, index):
    """
    Distinguish model item, given as link node information
    """
    if "src_vars()"==type or "dest_vars()"==type:
      return self._vars[index]
    if "src_cons()"==type:
      return self._cons_NL_all[index]
    if "src_objs()"==type:
      return self._objs_NL[index]
    assert not type.startswith("dest_cons(")      ## No grouping currently
    if "dest_objs()"==type:
      return self._objs[index]
    if type not in self._cons_Flat:
      raise Exception("Unknown constraint type: " + type)
    if index > len(self._cons_Flat[type]):
      raise Exception("Wrong index " + str(index) + " for constraint type " + type)
    return self._cons_Flat[type][index]


  # Match keyword to the original model
  def MatchOrigModel(self, keyw):
    result = {}
    result["NL Variables"] = self._matchRecords(self._vars, keyw, "is_from_nl")
    result["NL Defined Variables"] = self._matchRecords(self._dvars, keyw)
    result["NL Objectives"] = self._matchRecords(self._objs_NL, keyw)
#    result["NL Constraints"] \
#      = self._matchRecords(self._cons_NL.get("All"), keyw)
    result["NL Nonlinear Constraints"] \
      = self._matchRecords(self._cons_NL.get("Nonlinear"), keyw)
    result["NL Linear Constraints"] \
      = self._matchRecords(self._cons_NL.get("Linear"), keyw)
    result["NL Logical Constraints"] \
      = self._matchRecords(self._cons_NL.get("Logical"), keyw)
    result["NL SOS1 Constraints"] \
      = self._matchRecords(self._cons_NL.get("SOS1"), keyw)
    result["NL SOS2 Constraints"] \
      = self._matchRecords(self._cons_NL.get("SOS2"), keyw)
    return result

  # Match keyword to the final model
  def MatchFinalModel(self, keyw):
    result = {}
    result["Variables"] = self._matchRecords(self._vars, keyw)
    result["Objectives"] = self._matchRecords(self._objs, keyw)
    ctypes = FlatConTypes()
    for ct, cv in sorted(self._cons_Flat.items()):
      result["Constraints '" + \
        (ctypes[ct] if ct in ctypes else ct) + "'"] \
        = self._matchRecords(self._cons_Flat[ct], keyw)
    return result

  # Add records containing keyword
  # @return array of strings
  def _matchRecords(self, cnt, keyw, keyNeed1=None):
    result = ""
    if cnt is None:
      return result
    for i in cnt:
      if "final" not in i or 1==i["final"]:
        pr = str(i)         ## TODO printed form
        if "printed" in i:
          pr = i["printed"]
        assert len(pr)
        if ';'!=pr[-1]:
          pr = pr + ';'
        if (""==keyw or keyw in pr) \
            and (keyNeed1==None \
            or (keyNeed1 in i and 1==i[keyNeed1])):
          result = result + "  \n" + pr        ## Markdown: 2x spaces + EOL
    return result
