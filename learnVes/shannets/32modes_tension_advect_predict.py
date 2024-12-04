import torch
from Net_ves_advten_downsample_onelevel import Net_ves_advten_downsample_onelevel # from file import model class
import numpy as np


N = 32
modes = np.concatenate((np.arange(0, int(N/2)), np.arange(-int(N/2), 0)))
mod_list = np.where(np.abs(modes) <= 32)[0] + 1 # keeps indices in MATLAB order

nmodes = np.size(mod_list)# skip zeroth mode
model = Net_ves_advten_downsample_onelevel(14, 2.4, 30)

output_list = []
for ij in np.arange(0,nmodes-1):
  imode = mod_list[ij+1] # or [ij+1]
  imode_net = "/Users/gokberk/Documents/GitHub/ves2d/learnVes/shannets/32modes_advten_trained_models/2024Nov_downsample32_ves_advten_mode" + str(imode) + ".pth"
  model.load_state_dict(torch.load(imode_net, map_location="cpu"))
  model.eval()
  
  input_net = torch.from_numpy(input_shape[ij]).float()
  
  with torch.no_grad():
    output_net = model(input_net)
  output_list.append(output_net.detach().numpy())
  
     
  