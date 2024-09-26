import torch
from Ves_near_Fourier_downsample32 import Ves_near_Fourier_downsample32 # from file import model class
import numpy as np

# vesicle coordinates are normalized
# convert MATLAB's numpy array into PyTorch tensor
input_shape = torch.from_numpy(input_shape).float()

model = Ves_near_Fourier_downsample32(12, 1.6, 18)

output_list = []
N = 32
modes = np.concatenate((np.arange(0, int(N/2)), np.arange(-int(N/2), 0)))
mod_list = np.where(np.abs(modes) <= modesInUse)[0] + 1 # mod indices in MATLAB ordering

nmodes = np.size(mod_list)

for ij in np.arange(0,nmodes):
  imode = mod_list[ij]
  imode_net = "/Users/gokberk/Documents/GitHub/ves2d/learnVes/shannets/near_vel_32modesfft_models/Ves_downsample_nearFourier_nocoords_mode" + str(imode) + ".pth"
  model.load_state_dict(torch.load(imode_net, map_location="cpu"))
  model.eval()
  
  with torch.no_grad():
    output_net = model(input_shape)
  output_list.append(output_net.detach().numpy())
  
     
  