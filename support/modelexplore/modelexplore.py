
import streamlit as st
import json

from scripts.python.modelreader import ReadExplorerModel
from scripts.python.matcher import MatchSubmodel

# To work with local files in st 1.30.1, see
# https://discuss.streamlit.io/t/axioserror-request-failed-with-status-code-403/38112/13.
# The corresponding settings should not be used on a server.
uploader = st.sidebar.file_uploader(
    "Model file (JSONL)",
    help = "Reformulation file obtained by option `writegraph`  \n" + \
           "(https://mp.ampl.com/modeling-tools.html#reformulation-graph)")      # 2 spaces + EOL

# You can use a column just like st.sidebar:
srch = st.sidebar.text_input(
    'Search pattern:',
    help = "Pattern to filter the models' lines.\nLeave blank to see complete models.")

fwd = st.sidebar.checkbox(
    'Add descendants', disabled=True,
    help = 'When searching, include all solver model items derived from matching items in the NL model')
bwd = st.sidebar.checkbox(
    'Add ancestors', disabled=True,
    help = 'When searching, Include all NL/intermediate model items reduced to matching items in the solver model')

nl_ref_tree = st.sidebar.radio(
            "**NL model presentation mode:**",
            ["Text", "Reformulation tree :sparkles:"],
            help = "To best view the reformulation tree, download NL model as JSON and use a JSON viewer")
reftree = "Text"!=nl_ref_tree

left_column, right_column = st.columns(2)


# Cache the reading function
@st.cache_data
def ReadModel(uploader):
  return ReadExplorerModel(uploader)

# Cache the matching function?
# @st.cache_data  Need cacheable Model.
def MatchSelection(m, srch, fwd, bwd, reftree):
  return MatchSubmodel(m, srch, fwd, bwd, reftree)

# Write dictionary of entries
@st.cache_data
def WriteDict(d, title, reftree=False):
  whole = ("# "+title+"\n") if not reftree \
    else { "Title": title }         ## dict/json: add only non-empty sections
  for k, v in d.items():
    nv = v.count('\n') if not reftree else len(v)
    if nv:
      k1 = k  + ' (' + str(nv) + ')'
      if reftree:
        whole[k1] = v
      else:
        whole = whole + '\n\n##  ' + k1 + '\n'
        whole = whole + v
      with st.expander("""### """ + k1):
            with st.container(height=200):
              if reftree:
                st.json(v)
              else:
                st.code(v, language='ampl')
  return whole if not reftree else json.dumps(whole, indent=2)


filename_upl = ""
modelNL = ""
modelFlat = ""

# Or even better, call Streamlit functions inside a "with" block:
if uploader is not None:
  model = ReadModel(uploader)
  filename_upl = uploader.name
  subm1, subm2 = MatchSelection(model, srch, fwd, bwd, reftree)
  bytes1_data = subm1.GetData()
  bytes2_data = subm2.GetData()
  ttlEnd = " for '" + filename_upl + "' (search pattern: '" + srch + "')"
  with left_column:
    st.header("NL model",
        help = 'NL model lines matching the search pattern')
    st.write("Display mode:  **" + nl_ref_tree + "**")
    modelNL = WriteDict(bytes1_data, "NL Model"+ttlEnd, reftree)
  with right_column:
    st.header("Solver model",
        help = 'Solver model lines matching the search pattern')
    modelFlat = WriteDict(bytes2_data, "Solver model"+ttlEnd)
else:
  st.header("AMPL MP Reformulation Explorer")
  st.write("Documentation: https://mp.ampl.com/modeling-tools.html#reformulation-graph")
  st.divider()
  st.write("No file selected.")


st.sidebar.download_button("Download NL Model",
                   modelNL,
                   filename_upl + ('_NL.mod' if not reftree else '_NL.json'),
                   help = 'Download current NL model',
                   disabled = ("" == modelNL))
st.sidebar.download_button("Download Solver Model",
                   modelFlat,
                   filename_upl + '_solver.mod',
                   help = 'Download current solver model',
                   disabled = ("" == modelFlat))
