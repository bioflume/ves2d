import torch
from Net_ves_relax import Net_ves_relax # from file import model class

if __name__ == '__main__':

  model = Net_ves_relax(num_blocks=12)
  model.load_state_dict(torch.load("./Ves_relax.pth", map_location="cpu"))
#   model = torch.load('./ves_relax.pth')
  model.eval()
  dummy_input = torch.randn(1,2,128)
  
  traced_script_module = torch.jit.trace(model, (dummy_input))
  traced_script_module.save("ves_relax.pt")
  

#print(" ") 
#print('Model has been converted to PT') 