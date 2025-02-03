# -*- coding: utf-8 -*-

"""
Created on Fri Jan 13 23:57:34 2023

@author: sahiz
"""
import os
os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"
import torch
import matplotlib.pyplot as plt
import numpy as np
import time
from torchdiffeq import odeint,odeint_adjoint
from tqdm import tqdm
from PIL import Image
import glob
import os
import aux_fun as fn
if not os.path.exists("./plots"):
    os.makedirs("./plots")
#%% Data generation
"Load data"

dataSet = 'spiral'
tf = 2.5
dataSize = 200
y0 = [-5. , -5., -5.]

dataSet = 'lorentz'
tf = 15
dataSize = 400
y0 = [-5. , -5., 40.]
fn.generateData(y0=y0,dataSize=dataSize,name=dataSet,tf=tf)

data = np.float32(np.expand_dims(np.load('states_'+dataSet+'.npy').T,axis=1))
data_plot = data[:,0,:].T
t = np.load('t_'+dataSet+'.npy')
fn.visualize(t, data_plot)

#%% Model initialisation and training


"Define dynamic function"
n_hidden = 128

class ODEFunc(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.net = torch.nn.Sequential(torch.nn.Linear(3,n_hidden),
                                       # torch.nn.Tanh(),
                                       torch.nn.ReLU(),
                                       torch.nn.Linear(n_hidden,n_hidden),
                                       torch.nn.ReLU(),
                                       torch.nn.Linear(n_hidden,n_hidden),
                                       torch.nn.ReLU(),
                                       torch.nn.Linear(n_hidden,3))
        for m in self.net.modules():
            if isinstance(m,torch.nn.Linear):
                torch.nn.init.normal_(m.weight, mean = 0 , std = 0.1)
                torch.nn.init.constant_(m.bias, val=0)
    
    def forward(self, t, y):
        output = self.net(y)
        return output
    
def getBatch(batchSize, batchLength):
    start = torch.from_numpy(np.random.choice(np.arange(dataSize - batchLength, dtype=np.int64), batchSize, replace = False))
    initial_states = observations[start]
    times  = t_torch[:batchLength]
    targets  = torch.stack([observations[start + i] for i in range(batchLength)],dim = 0)
    return initial_states, times, targets

    
"Specify observations and batch length"
observations = torch.from_numpy(data)
t_torch = torch.from_numpy(t)
initial_state = observations[0]
dataSize = len(t)
batchLength = 15
batchSize = dataSize//batchLength

"Initialize model and parameters"
func = ODEFunc()
loss_fn = torch.nn.MSELoss()
optimizer = torch.optim.Adam(func.parameters(), lr=1e-2)

epochs = 2000
losses = []
losses_val = []

startTime = time.time()

for epoch in tqdm(range(epochs+1)):
    if epoch > 800:
        optimizer.param_groups[0]['lr'] = 1e-3
    optimizer.zero_grad()
    
    initial_states, times, targets = getBatch(batchSize, batchLength)
    
    pred_Y = odeint_adjoint(func = func, y0=initial_states, t=times, rtol = 1e-7, atol =1e-9, method='dopri5')
    # pred_Y = odeint_adjoint(func = func, y0=initial_states, t=times, rtol = 1e-2, atol =1e-3, method='dopri5')
    loss = loss_fn(pred_Y,targets)
    loss.backward()
    optimizer.step()
    losses.append(loss.item())
    if epoch % 50 == 0:
        with torch.no_grad():
            prediction = odeint(func,initial_state, t_torch, rtol = 1e-7, atol =1e-9, method='dopri5')
            loss = loss_fn(prediction, observations)
            prediction = np.array(prediction).T.squeeze()
            fn.visualize(t,data_plot,prediction,iter=epoch,save='plots/NODE '+dataSet+ f' prediction epoch {epoch}')
            print(f'Iteration {epoch} | Total loss {loss.item():.2e}')
            losses_val.append(loss.item())
endTime = time.time() - startTime

print('----------------------------')
print(f'Training done in {endTime/60:.1f} min')
plt.plot(losses,label='training loss')
plt.plot(np.arange(epochs,step=49),losses_val,label='validation loss')
plt.yscale('log')
plt.grid()
plt.xlabel('epochs')
plt.ylabel('error')
plt.title('error evolution')
plt.legend()
plt.tight_layout()
plt.show()

fn.make_gif(name='NODE '+dataSet)

"Save the model"     
torch.save(func.state_dict(), 'NODE_model_state_dict_'+dataSet)


#%% Test prediction on unseen data

if dataSize<250:
    fn.generateData(y0=data[-1,0,:],dataSize=3*dataSize,name=dataSet,t0=t[-1],tf=t[-1]+3*tf)
else:
    fn.generateData(y0=data[-1,0,:],dataSize=dataSize,name=dataSet,t0=t[-1],tf=t[-1]+tf)
data_forecast = np.float32(np.expand_dims(np.load('states_'+dataSet+'.npy').T,axis=1))
data_plot = data_forecast[:,0,:].T

initial_state_forecast = observations[-1,:,:]
forecast = torch.from_numpy(data_forecast)
t_forecast = torch.from_numpy(np.load('t_'+dataSet+'.npy'))


with torch.no_grad():
    prediction = odeint(func,initial_state_forecast, t_forecast, rtol = 1e-7, atol =1e-9, method='dopri5')
    loss = loss_fn(prediction, forecast)
    print(f'Forecast error is {loss.item():.2e}')
prediction = np.array(prediction).T.squeeze()
fn.visualize(t_forecast,data_plot,prediction,iter='Forecast',save='NODE forecast '+dataSet)