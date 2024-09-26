import torch
from Ves_relax_downsample_zerolevel import  Ves_relax_downsample_zerolevel # from file import model class

# convert MATLAB's numpy array into PyTorch tensor
input_shape = torch.from_numpy(input_shape).float()

model = Ves_relax_downsample_zerolevel(12, 2.5)
model.load_state_dict(torch.load("/Users/gokberk/Documents/GitHub/ves2d/learnVes/shannets/Ves_relax_1e-5_downsample_DIFF.pth", map_location="cpu"))
model.eval()
  
predicted_shape = (model(input_shape))
predicted_shape = predicted_shape.detach().numpy()
     