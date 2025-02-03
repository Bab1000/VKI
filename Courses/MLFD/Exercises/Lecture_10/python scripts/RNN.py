# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 21:56:56 2023

@author: sahiz
"""
import os
os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"

import torch
from torch import nn
import numpy as np
import matplotlib.pyplot as plt
import aux_fun as fn
from time import time

if not os.path.exists("./plots"):
    os.makedirs("./plots")
    
#%% Data Generation
"Load data"

dataSet = 'spiral'
tf = 2.5
dataSize = 200
y0 = [-5. , -5., -5.]

# dataSet = 'lorentz'
# tf = 15
# dataSize = 400
# y0 = [-5. , -5., 40.]

fn.generateData(y0=y0,dataSize=dataSize,name=dataSet,tf=tf)

data = np.float32(np.expand_dims(np.load('states_'+dataSet+'.npy').T,axis=1))
data_plot = data[:,0,:].T
t = np.load('t_'+dataSet+'.npy')
fn.visualize(t, data_plot)
#%% Model init and training
class RNN(nn.Module):
  def __init__(self, input_size, hidden_dim, output_size, n_layers):
    super(RNN, self).__init__()
    self.hidden_dim = hidden_dim
    self.n_layers = n_layers
    self.rnn = nn.RNN(input_size, 
                      hidden_dim, n_layers, 
                      # nonlinearity='relu',
                      nonlinearity='tanh',
                      batch_first=True) # RNN hidden units
    self.fc = nn.Linear(hidden_dim, output_size) # output layer

  def forward(self, x):
    bs, _, _ = x.shape
    h0 = torch.zeros(self.n_layers, bs, self.hidden_dim).requires_grad_()
    h_n, hidden = self.rnn(x, h0.detach())
    h_n = h_n.view(bs, -1, self.hidden_dim)
    out = self.fc(h_n)
    return out[:, -1, :]




def create_batch(input_data, batchLength):
    train_data = []
    target = []
    L = len(input_data)
    for i in range(L-batchLength):
        sequence = input_data[i:i+batchLength]
        nextStep = input_data[i+batchLength:i+batchLength+1]
        train_data.append(torch.FloatTensor(sequence).squeeze(1).unsqueeze(0))
        target.append(torch.FloatTensor(nextStep).squeeze(1).unsqueeze(0))

    return torch.cat(train_data, 0), torch.cat(target, 0)

"Split data into training and targets"
batchLength = 15
observations = torch.from_numpy(data)
data_train, targets = create_batch(observations, batchLength)
data_train.shape
epochs = 2000
losses = []
losses_val = []

"Initialize model and parameters"
model = RNN(input_size=3, hidden_dim=128, output_size=3, n_layers=1)
loss_function = nn.MSELoss()
optimizer = torch.optim.Adam(model.parameters(), lr=1e-2)

"training loop"
begin = time()
for epoch in range(epochs+1):
    if epoch > 800:
        optimizer.param_groups[0]['lr'] = 1e-3
    "Training"
    optimizer.zero_grad()
    y_pred = model(data_train)
    loss = loss_function(y_pred, targets.squeeze())
    loss.backward()
    optimizer.step()
    losses.append(loss.item())
    
    "Validation"
    if epoch%50 == 0:
        prediction = 0
        test_inputs = data_train[0].detach().numpy()
        prediction = []
        for i in range(batchLength):
            prediction.append(test_inputs[i])

        with torch.no_grad():
            for i in range(len(t)-batchLength):
                seq = torch.tensor(np.array(prediction).squeeze()[-batchLength:])
                model_out = model(seq.unsqueeze(0)).squeeze().detach().numpy().reshape(-1)
                prediction.append(model_out)
        prediction = np.array(prediction).T
        fn.visualize(t,data_plot,prediction,iter=epoch,save='plots/RNN '+dataSet+ f' prediction epoch {epoch}')
        prediction = torch.from_numpy(prediction)
        loss = loss_function(prediction.T,observations[:,0,:]).item()
        losses_val.append(loss)
        print(f'Epoch {epoch}')
        print(f'Loss on the complete trajectory: {loss:.2e}')
        
end = time() - begin
print('----------------------------')
print(f'Training done in {end/60:.1f} min')
plt.plot(losses,label='training loss')
plt.plot(np.arange(epochs+50,step=50),losses_val,label='validation loss')
plt.yscale('log')
plt.grid()
plt.xlabel('epochs')
plt.ylabel('error')
plt.title('error evolution')
plt.legend()
plt.tight_layout()
plt.show()

fn.make_gif(name='RNN '+dataSet)


#%% Assess forecasting capabilities
if dataSize<250:
    fn.generateData(y0=data[-1,0,:],dataSize=3*dataSize,name=dataSet,t0=t[-1],tf=t[-1]+3*tf)
else:
    fn.generateData(y0=data[-1,0,:],dataSize=dataSize,name=dataSet,t0=t[-1],tf=t[-1]+tf)
data = np.float32(np.expand_dims(np.load('states_'+dataSet+'.npy').T,axis=1))
data_plot = data[:,0,:].T
t = np.load('t_'+dataSet+'.npy')
fn.visualize(t, data_plot)

forecast_inputs = observations.squeeze().detach().numpy()
prediction = []
for i in range(batchLength):
    prediction.append(forecast_inputs[-batchLength+i])

with torch.no_grad():
    for i in range(len(t)):
        seq = torch.tensor(np.array(prediction).squeeze()[-batchLength:])
        model_out = model(seq.unsqueeze(0)).squeeze().detach().numpy().reshape(-1)
        prediction.append(model_out)
prediction = np.array(prediction).T
fn.visualize(t,data_plot,prediction[:,batchLength:],iter= 'Forecast',save='Forecast RNN '+dataSet)
loss = loss_function(torch.from_numpy(prediction[:,batchLength:]),torch.from_numpy(data_plot))
print(f'Forecast error is {loss:.2e}')

