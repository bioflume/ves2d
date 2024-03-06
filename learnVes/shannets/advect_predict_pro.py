import torch
from Net_ves_adv_fft import Net_ves_adv_fft # from file import model class
import numpy as np

# vesicle coordinates are normalized
# convert MATLAB's numpy array into PyTorch tensor
input_shape = torch.from_numpy(input_shape).float()

theta = np.arange(128)/128*2*np.pi
theta = theta.reshape(128,1)
bases = 1/128*np.exp(1j*theta*np.arange(128).reshape(1,128))  

rrMat = np.real(bases)
iiMat = np.imag(bases)
  
rr_ii_Mat = np.concatenate((rrMat, iiMat)).transpose() 
# [256, 128] shape transposed to 128 x 256
rr_ii_Mat = torch.from_numpy(rr_ii_Mat).float()


model = Net_ves_adv_fft(12, 1.7, 20)

output_list = []
for imode in np.arange(2,65):
  imode_net = "/Users/gokberkkabacaoglu/Documents/GitHub/ves2d/learnVes/shannets/ves_fft_models/ves_fft_mode" + str(imode) + ".pth"
  model.load_state_dict(torch.load(imode_net, map_location="cpu"))
  model.eval()

  # Get the fourier mode corresponding to imode
  basis = rr_ii_Mat[imode-1,:].reshape(1,1,256).float()
    
  # tile the fourier mode array for num_ves
  fourX_tile = torch.tile(basis,(num_ves,1,1))
  
  input_net = torch.zeros(num_ves,2,256)
  input_net[:,0,:128] = input_shape[:,0,:]
  input_net[:,0,128:] = input_shape[:,1,:]
  input_net[:,1,:] = fourX_tile
  
  with torch.no_grad():
    output_net = model(input_net)
  output_list.append(output_net.detach().numpy())
  
     
  