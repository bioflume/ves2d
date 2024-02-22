import torch
import numpy as np
from Net_ves_relax_midfat import Net_ves_midfat # from file import model class

if __name__ == '__main__':
  num_ves = 1
  input_shape = torch.randn(1,2,128)
  model = Net_ves_midfat(num_blocks=14)
  model.load_state_dict(torch.load("./ves_relax160k_dt1e-5.pth", map_location="cpu"))
  model.eval()
  # input shape is (N,2,128) N is the number of vesicles input
  
  # normalize the input
  # input normalizing parameters
  x_mean = 2.8696319626098088e-11 
  x_std = 0.06018252670764923
  y_mean = 2.8696319626098088e-11 
  y_std = 0.13689695298671722
  

  for ik in range(num_ves):
    input_shape[ik,0,:] = (input_shape[ik,0,:] - x_mean)/x_std
    input_shape[ik,1,:] = (input_shape[ik,1,:] - y_mean)/y_std
  
  predicted_shape = model(input_shape)
  
  # normalize output
  x_mean = 1.8570537577033974e-05
  x_std = 0.058794111013412476
  y_mean = -0.0005655255517922342 
  out_y_std = 0.13813920319080353
 
  for ik in range(num_ves):
    predicted_shape[ik,0,:] = x_std*predicted_shape[ik,0,:] + x_mean
    predicted_shape[ik,1,:] = y_std*predicted_shape[ik,1,:] + y_mean
  
  print(predicted_shape)
  
  