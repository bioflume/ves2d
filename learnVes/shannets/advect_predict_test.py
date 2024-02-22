import torch
import numpy as np
from Net_ves_adv_fft import Net_ves_adv_fft # from file import model class

if __name__ == '__main__':
  in_param = np.load('ves_fft_in_param.npy')
  out_param = np.load('ves_fft_out_param.npy')
  
  num_ves = 1
  input_shape = torch.randn(1,2,128)
  model = Net_ves_adv_fft(12, 1.7, 20)
  # fourier mode
  theta = np.arange(128)/128*2*np.pi
  theta = theta.reshape(128,1)
  bases = 1/128*np.exp(1j*theta*np.arange(128).reshape(1,128))  
  
  output_list = []
  for imode in np.arange(2,129):

    imode_net = "./ves_fft_models/ves_fft_mode" + str(imode) + ".pth"
    model.load_state_dict(torch.load(imode_net, map_location="cpu"))
    model.eval()
    # input shape is (N,2,256) N is the number of vesicles input
    # input normalizing parameters
    x_mean = in_param[imode-2,0]
    x_std = in_param[imode-2,1]
    y_mean = in_param[imode-2,2]
    y_std = in_param[imode-2,3]
    
    # Get the fourier mode corresponding to imode
    rr, ii = np.real(bases[:,imode-1]), np.imag(bases[:,imode-1])
    basis = torch.from_numpy(np.concatenate((rr,ii))).float().reshape(1,1,256).float()
    
    # tile the fourier mode array for num_ves
    fourX_tile = torch.tile(basis,(num_ves,1,1))
    
    # normalize input and add the associated fourier mode vectors
    input_net = torch.zeros(num_ves,2,256)
    input_net[:,0,:128] = (input_shape[:,0,:] - x_mean)/x_std
    input_net[:,0,128:] = (input_shape[:,1,:] - y_mean)/y_std
    input_net[:,1,:] = fourX_tile
    
    # predict
    output_net = model(input_net)
  
    # denormalize output
    real_mean = out_param[imode-2,0]
    real_std = out_param[imode-2,1]
    imag_mean = out_param[imode-2,2]
    imag_std = out_param[imode-2,3]
    
    output_net[:,0,:] = real_std*output_net[:,0,:] + real_mean
    output_net[:,1,:] = imag_std*output_net[:,1,:] + imag_mean
    
    output_list.append(output_net.detach().numpy())
    print(output_net)
    input()
  
  