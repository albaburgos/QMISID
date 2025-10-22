import torch
import torch.nn as nn
import torch.nn.functional as F
import matplotlib.pyplot as plt

class KappaNet(nn.Module):
    def __init__(self, input_dim, hidden_dim=64):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, 1)
        )

    def forward(self, e):  
        return self.net(e).squeeze(-1)  

def f_phi(phi1, phi2):
    term1 = F.softplus(phi1 + phi2)
    stacked = torch.stack([phi1, phi2], dim=-1)  
    term2 = torch.logsumexp(stacked, dim=-1)    
    return term1 - term2

def pairwise_loss(model, e1, e2, labels):
    phi1, phi2 = model(e1), model(e2)           
    f = f_phi(phi1, phi2)                       
    return F.binary_cross_entropy_with_logits(f, labels.float())

device = "cpu" 
D, N = 10, 32
model = KappaNet(input_dim=D).to(device)
e1 = torch.randn(N, D, device=device)
e2 = torch.randn(N, D, device=device)
labels = torch.randint(0, 2, (N,), device=device)

optimizer = torch.optim.Adam(model.parameters(), lr= 1e-3)

for epoch in range(1000):
    loss = pairwise_loss(model, e1, e2, labels)
    optimizer.zero_grad()
    loss.backward()
    optimizer.step()
    print("loss:", float(loss))


losses = []

for epoch in range(1000):
    loss = pairwise_loss(model, e1, e2, labels)
    optimizer.zero_grad()
    loss.backward()
    optimizer.step()
    losses.append(loss.item())

plt.figure(figsize=(8, 5))
plt.plot(losses, label='Training Loss')
plt.xlabel('Epoch')
plt.ylabel('Loss')
plt.title('KappaNet Training Loss Over Time')
plt.legend()
plt.grid(True)
plt.show()