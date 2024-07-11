import torch
from Net_ves_advten import Net_ves_advten # from file import model class
import numpy as np


N = 128
modes = np.concatenate((np.arange(0, int(N/2)), np.arange(-int(N/2), 0)))
mod_list = np.where(np.abs(modes) <= 16)[0] + 1 # keeps indices in MATLAB order

nmodes = np.size(mod_list)# skip zeroth mode
model = Net_ves_advten(12, 2.5, 24)

output_list = []
for ij in np.arange(0,nmodes-1):
  imode = mod_list[ij] # or [ij+1]
  imode_net = "/Users/gokberk/Documents/GitHub/ves2d/learnVes/shannets/ves_advten_models/ves_advten_mode" + str(imode) + ".pth"
  model.load_state_dict(torch.load(imode_net, map_location="cpu"))
  model.eval()
  
  input_net = torch.from_numpy(input_shape[ij]).float()
  
  with torch.no_grad():
    output_net = model(input_net)
  output_list.append(output_net.detach().numpy())
  
     
  