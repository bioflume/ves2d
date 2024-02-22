import torch
from Net_ves_relax import Net_ves_relax # from file import model class

# convert MATLAB's numpy array into PyTorch tensor
input_shape = torch.from_numpy(input_shape).float()

model = Net_ves_relax(num_blocks=12)
model.load_state_dict(torch.load("/Users/gokberk/Documents/GitHub/ves2d/learnVes/shannets/ves_relax.pth", map_location="cpu"))
model.eval()
with torch.no_grad():
    predicted_shape = (model(input_shape))
    
predicted_shape = predicted_shape.detach().numpy()
     
  