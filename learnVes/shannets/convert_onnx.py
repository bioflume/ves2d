import torch
from Net_ves_relax import Net_ves_relax # from file import model class

if __name__ == '__main__':

  model = Net_ves_relax(num_blocks=12)
  model.load_state_dict(torch.load("./Ves_relax.pth", map_location="cpu"))
#   model = torch.load('./ves_relax.pth')
  model.eval()
  dummy_input = torch.randn(1,2,128)
  input_names = ["actual_input"]
  output_names = ["output"]

  torch.onnx.export(model, 
                    dummy_input,
                    "ves_relax.onnx",
                    verbose=False,
                    input_names=input_names,
                    output_names=output_names,
                    export_params=True,
                    dynamic_axes={'actual_input' : {0 : 'batch_size'},    # variable length axes 
                                'output' : {0 : 'batch_size'}}
                    )
print(" ") 
print('Model has been converted to ONNX') 