import torch
from Net_ves_adv_fft import Net_ves_adv_fft # from file import model class

#import numpy as np
# Input all the vesicles
# vesicle coordinates are normalized
# convert MATLAB's numpy array into PyTorch tensor
input_shape = torch.from_numpy(input_shape).float()
# input_shape contains num_ves 
#theta = np.arange(128)/128*2*np.pi
#theta = theta.reshape(128,1)
#bases = 1/128*np.exp(1j*theta*np.arange(128).reshape(1,128))  
  
model = Net_ves_adv_fft(12, 1.7, 20)

output_list = []
for imode in modeList:
  imode_net = "./ves_fft_models/ves_fft_mode" + str(imode) + ".pth"
  model.load_state_dict(torch.load(imode_net, map_location="cpu"))
  model.eval()
    
  # Get the fourier mode corresponding to imode
  #rr, ii = np.real(bases[:,imode-1]), np.imag(bases[:,imode-1])
  #basis = torch.from_numpy(np.concatenate((rr,ii))).float().reshape(1,1,256).float()
  basis = torch.from_numpy(basesConcat[:,imode-1]).float().reshape(1,1,256).float()
  # tile the fourier mode array for num_ves
  fourX_tile = torch.tile(basis,(num_ves,1,1))
  
  input_net = torch.zeros(num_ves,2,256)
  input_net[:,0,:128] = input_shape[:,0,:]
  input_net[:,0,128:] = input_shape[:,1,:]
  input_net[:,1,:] = fourX_tile
  
  output_net = (model(input_net))
  output_list.append(output_net.detach().numpy())
  
     
  