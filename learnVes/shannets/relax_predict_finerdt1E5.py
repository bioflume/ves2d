import torch
from Net_ves_relax_midfat import Net_ves_midfat # from file import model class

# convert MATLAB's numpy array into PyTorch tensor
input_shape = torch.from_numpy(input_shape).float()

model = Net_ves_midfat(num_blocks=16)
model.load_state_dict(torch.load("/Users/gokberk/Documents/GitHub/ves2d/learnVes/shannets/ves_relax619k_mirr_finer_dt1e-5.pth", map_location="cpu"))
model.eval()
  
predicted_shape = (model(input_shape))
predicted_shape = predicted_shape.detach().numpy()
     