import torch
import torch.nn as nn
import math

# input: coords 2ch (2,128)
# starter: in 2ch, out 12ch= 10ch + coords 2ch, coords is duplicated
# evo: in 12ch, out 12ch = 10ch + coords 2ch
# end: in 12ch, out 6ch 

class Starter_Block(nn.Module):
    def __init__(self, factor, ch):
        super(Starter_Block, self).__init__()
        conv1_out_ch = math.floor(9*factor)
        self.conv1_1 = nn.Conv1d(2, conv1_out_ch, kernel_size=7, stride=1, padding=3, padding_mode='circular')
        self.conv1_2 = nn.Conv1d(2, conv1_out_ch, kernel_size=7, stride=1, padding=6, dilation=2, padding_mode='circular')
        self.bn1 = nn.BatchNorm1d(num_features=2*conv1_out_ch)
        self.relu = nn.ReLU(inplace=True)
        conv2_out_ch = math.floor(10*factor)
        self.conv2_1 = nn.Conv1d(2*conv1_out_ch, conv2_out_ch, kernel_size=7, stride=1, padding=3, padding_mode='circular')
        self.conv2_2 = nn.Conv1d(2*conv1_out_ch, conv2_out_ch, kernel_size=7, stride=1, padding=6, dilation=2, padding_mode='circular')
        deconv_ch = 2*conv2_out_ch
        self.bn2 = nn.BatchNorm1d(num_features=deconv_ch)
        # self.deconv1 = nn.ConvTranspose1d(deconv_ch, deconv_ch, kernel_size=10, stride=1, padding=4)
        # self.bn3 = nn.BatchNorm1d(num_features=deconv_ch)
        # self.deconv2 = nn.ConvTranspose1d(deconv_ch, deconv_ch, kernel_size=10, stride=1, padding=4)
        # self.bn4 = nn.BatchNorm1d(num_features=deconv_ch)
        self.conv3_1 = nn.Conv1d(2*conv2_out_ch, conv2_out_ch, kernel_size=7, stride=1, padding=3, padding_mode='circular')
        self.conv3_2 = nn.Conv1d(2*conv2_out_ch, conv2_out_ch, kernel_size=7, stride=1, padding=6, dilation=2, padding_mode='circular')
        self.bn3 = nn.BatchNorm1d(num_features=deconv_ch)

        self.conv4_1 = nn.Conv1d(2*conv2_out_ch, conv2_out_ch, kernel_size=7, stride=1, padding=3, padding_mode='circular')
        self.conv4_2 = nn.Conv1d(2*conv2_out_ch, conv2_out_ch, kernel_size=7, stride=1, padding=6, dilation=2, padding_mode='circular')
        self.bn4 = nn.BatchNorm1d(num_features=deconv_ch)

        self.mix = nn.Conv1d(deconv_ch, ch, kernel_size=3, stride=1, padding=1)
        # self.W = nn.Conv1d(1, 1, kernel_size=3, stride=1, padding=1)

    def forward(self, input):

        out_1_1 = self.relu(self.conv1_1(input))
        out_1_2 = self.relu(self.conv1_2(input))
        out = torch.concat((out_1_1, out_1_2), dim=1)
        out = self.bn1(out)

        out_2_1 = self.relu(self.conv2_1(out))
        out_2_2 = self.relu(self.conv2_2(out))
        out = torch.concat((out_2_1, out_2_2), dim=1)
        out = self.bn2(out)

        # out = self.relu(self.deconv1(out))
        # out = self.bn3(out)
        # out = self.relu(self.deconv2(out))
        # out = self.bn4(out)

        out_3_1 = self.relu(self.conv3_1(out))
        out_3_2 = self.relu(self.conv3_2(out))
        out = torch.concat((out_3_1, out_3_2), dim=1)
        out = self.bn3(out)

        out_4_1 = self.relu(self.conv4_1(out))
        out_4_2 = self.relu(self.conv4_2(out))
        out = torch.concat((out_4_1, out_4_2), dim=1)
        out = self.bn4(out)

        out = self.mix(out)
        # out = torch.cat((out, input), dim=1)
        
        return out

class Evo_Block(nn.Module):
    def __init__(self, factor, ch):
        super(Evo_Block, self).__init__()
        conv1_out_ch = math.floor(8*factor)
        self.conv1_1 = nn.Conv1d(ch, conv1_out_ch, kernel_size=7, stride=1, padding=3, padding_mode='circular')
        self.conv1_2 = nn.Conv1d(ch, conv1_out_ch, kernel_size=7, stride=1, padding=6, dilation=2, padding_mode='circular')
        self.bn1 = nn.BatchNorm1d(num_features=2*conv1_out_ch)
        self.relu = nn.ReLU(inplace=True)
        conv2_out_ch = math.floor(10*factor)
        self.conv2_1 = nn.Conv1d(2*conv1_out_ch, conv2_out_ch, 
                                 kernel_size=7, stride=1, padding=3, padding_mode='circular')
        self.conv2_2 = nn.Conv1d(2*conv1_out_ch, conv2_out_ch, 
                                 kernel_size=7, stride=1, padding=6, dilation=2, padding_mode='circular')
        deconv_ch = 2*conv2_out_ch
        self.bn2 = nn.BatchNorm1d(num_features=deconv_ch)
        # self.deconv1 = nn.ConvTranspose1d(deconv_ch, deconv_ch, kernel_size=10, stride=1, padding=4)
        # self.bn3 = nn.BatchNorm1d(num_features=deconv_ch)
        # self.deconv2 = nn.ConvTranspose1d(deconv_ch, deconv_ch, kernel_size=10, stride=1, padding=4)
        # self.bn4 = nn.BatchNorm1d(num_features=deconv_ch)

        self.conv3_1 = nn.Conv1d(2*conv2_out_ch, conv2_out_ch, kernel_size=7, stride=1, padding=3, padding_mode='circular')
        self.conv3_2 = nn.Conv1d(2*conv2_out_ch, conv2_out_ch, kernel_size=7, stride=1, padding=6, dilation=2, padding_mode='circular')
        self.bn3 = nn.BatchNorm1d(num_features=deconv_ch)

        self.conv4_1 = nn.Conv1d(2*conv2_out_ch, conv2_out_ch, kernel_size=7, stride=1, padding=3, padding_mode='circular')
        self.conv4_2 = nn.Conv1d(2*conv2_out_ch, conv2_out_ch, kernel_size=7, stride=1, padding=6, dilation=2, padding_mode='circular')
        self.bn4 = nn.BatchNorm1d(num_features=deconv_ch)

        self.mix = nn.Conv1d(deconv_ch, ch, kernel_size=3, stride=1, padding=1)
        # self.W = nn.Conv1d(1, 1, kernel_size=3, stride=1, padding=1)

    def forward(self, input):

        out_1_1 = self.relu(self.conv1_1(input))
        out_1_2 = self.relu(self.conv1_2(input))
        out = torch.concat((out_1_1, out_1_2), dim=1)
 
        out = self.bn1(out)

        out_2_1 = self.relu(self.conv2_1(out))
        out_2_2 = self.relu(self.conv2_2(out))
        out = torch.concat((out_2_1, out_2_2), dim=1)
 
        out = self.bn2(out)

        out_3_1 = self.relu(self.conv3_1(out))
        out_3_2 = self.relu(self.conv3_2(out))
        out = torch.concat((out_3_1, out_3_2), dim=1)
        out = self.bn3(out)

        out_4_1 = self.relu(self.conv4_1(out))
        out_4_2 = self.relu(self.conv4_2(out))
        out = torch.concat((out_4_1, out_4_2), dim=1)
        out = self.bn4(out)

        out = self.mix(out)

        out = out + input
       
        return out

class End_Block(nn.Module):
    def __init__(self, factor, ch):
        super(End_Block, self).__init__()
        conv1_out_ch = math.floor(8*factor)
        self.conv1_1 = nn.Conv1d(ch, conv1_out_ch, kernel_size=7, stride=1, padding=3, padding_mode='circular')
        self.conv1_2 = nn.Conv1d(ch, conv1_out_ch, kernel_size=7, stride=1, padding=6, dilation=2, padding_mode='circular')
        self.bn1 = nn.BatchNorm1d(num_features=2*conv1_out_ch)
        self.relu = nn.ReLU(inplace=True)
        conv2_out_ch = math.floor(10*factor)
        self.conv2_1 = nn.Conv1d(2*conv1_out_ch, conv2_out_ch, 
                                 kernel_size=7, stride=1, padding=3, padding_mode='circular')
        self.conv2_2 = nn.Conv1d(2*conv1_out_ch, conv2_out_ch, 
                                 kernel_size=7, stride=1, padding=6, dilation=2, padding_mode='circular')
        deconv_ch = 2*conv2_out_ch
        self.bn2 = nn.BatchNorm1d(num_features=deconv_ch)
        # self.deconv1 = nn.ConvTranspose1d(deconv_ch, deconv_ch, kernel_size=10, stride=1, padding=4)
        # self.bn3 = nn.BatchNorm1d(num_features=deconv_ch)
        # self.deconv2 = nn.ConvTranspose1d(deconv_ch, deconv_ch, kernel_size=10, stride=1, padding=4)
        # self.bn4 = nn.BatchNorm1d(num_features=deconv_ch)
        # self.deconv3 = nn.ConvTranspose1d(deconv_ch, deconv_ch, kernel_size=10, stride=1, padding=4)
        # self.bn5 = nn.BatchNorm1d(num_features=deconv_ch)

        self.conv3_1 = nn.Conv1d(2*conv2_out_ch, conv2_out_ch, kernel_size=7, stride=1, padding=3, padding_mode='circular')
        self.conv3_2 = nn.Conv1d(2*conv2_out_ch, conv2_out_ch, kernel_size=7, stride=1, padding=6, dilation=2, padding_mode='circular')
        self.bn3 = nn.BatchNorm1d(num_features=deconv_ch)

        self.conv4_1 = nn.Conv1d(2*conv2_out_ch, conv2_out_ch, kernel_size=7, stride=1, padding=3, padding_mode='circular')
        self.conv4_2 = nn.Conv1d(2*conv2_out_ch, conv2_out_ch, kernel_size=7, stride=1, padding=6, dilation=2, padding_mode='circular')
        self.bn4 = nn.BatchNorm1d(num_features=deconv_ch)

        self.mix = nn.Conv1d(deconv_ch, 12, kernel_size=3, stride=1, padding=1)
        # self.W = nn.Conv1d(1, 1, kernel_size=3, stride=1, padding=1)

    def forward(self, input):
        # residual = input[:,:-2]
        # mode = input[:,-2:]

        out_1_1 = self.relu(self.conv1_1(input))
        out_1_2 = self.relu(self.conv1_2(input))
        out = torch.concat((out_1_1, out_1_2), dim=1)
        # del out_1_1
        # del out_1_2
        out = self.bn1(out)

        out_2_1 = self.relu(self.conv2_1(out))
        out_2_2 = self.relu(self.conv2_2(out))
        out = torch.concat((out_2_1, out_2_2), dim=1)
        # del out_2_1
        # del out_2_2
        out = self.bn2(out)


        out_3_1 = self.relu(self.conv3_1(out))
        out_3_2 = self.relu(self.conv3_2(out))
        out = torch.concat((out_3_1, out_3_2), dim=1)
        out = self.bn3(out)

        out_4_1 = self.relu(self.conv4_1(out))
        out_4_2 = self.relu(self.conv4_2(out))
        out = torch.concat((out_4_1, out_4_2), dim=1)
        out = self.bn4(out)

        out = self.mix(out)

        # out_u = torch.concat((out, mode), dim=1)
       
        return out


class Ves_near_Fourier_downsample32(nn.Module):
    def __init__(self, num_blocks, factor, ch):
        super(Ves_near_Fourier_downsample32, self).__init__()
        self.factor = factor
        self.ch = ch
        self.layer = self.make_layer(num_blocks-1, ch)
        self.starter = Starter_Block(factor, ch)
        self.end = End_Block(factor, ch)
        
    def make_layer(self, num_blocks, ch):
        layers = []
        for _ in range(num_blocks):
            layers.append(Evo_Block(self.factor, ch))
        return nn.Sequential(*layers)
    
    def forward(self, input):
        out = self.layer(self.starter(input))
        ans = self.end(out)
        return ans








