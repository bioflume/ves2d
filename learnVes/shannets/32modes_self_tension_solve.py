import torch
from Ves_selften_downsample_zerolevel import pdeNet_Ves_fat_factor_modein6_zerolevel # from file import model class

# convert MATLAB's numpy array into PyTorch tensor
input_shape = torch.from_numpy(input_shape).float()

model = pdeNet_Ves_fat_factor_modein6_zerolevel(12, 1.5, 20)
model.load_state_dict(torch.load("/Users/gokberk/Documents/GitHub/ves2d/learnVes/shannets/ves_downsample_selften_zerolevel.pth", map_location="cpu"))
model.eval()
  
predicted_shape = (model(input_shape))
predicted_shape = predicted_shape.detach().numpy()
     