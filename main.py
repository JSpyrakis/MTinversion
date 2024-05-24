import numpy as np
import pandas as pd
from statsmodels.regression.linear_model import OLS
from statsmodels.formula.api import mixedlm
from statsmodels.tools.tools import add_constant
import statsmodels.api as sm 
import matplotlib.pyplot as plt
import netCDF4 as nc
import scipy as sc
import math 

# Allowing flexiblity in variables
timepoints = 180
NumStations = 6
NumStationComponents = [3,3,3,3,3,3] 
TotalComponents = sum(NumStationComponents)
N_total = 18*180

# When creating parameters in a more general case, this array may be an input.

DATA1 = pd.read_csv('merge1.csv')
SplitData = []
# Splitting DATA into mini arrays for each component
for i in range(sum(NumStationComponents)):
    SplitData.append(DATA1.iloc[i*180:(i+1)*180,:])

d = np.array([SplitData])



DATA = pd.read_csv('merge1.csv')
X = DATA[['g1', 'g2', 'g3', 'g4', 'g5', 'g6']]
DATA['value'] = DATA['value']
x_matrix = X.to_numpy()
#print(x_matrix)


# Conducting a Crude Least Squares for parameter estimation
Crude_LS = OLS(DATA['value'], X ).fit()
mm0 = Crude_LS
beta = Crude_LS.params.values
LS_RSE = (Crude_LS.resid.std(ddof=X.shape[1]))**2 #Finding the Residual Std Error Squared

# Constructing S matrix for timeshifts co domain
n = 5 # Adding variable option for later
s = np.arange(-n,n+1,1)

# Defining Parameters
m_group = 18
tau = 0
LL= 3
KL = 180/LL
N = N_total*len(s)
x_comb = []
for i in range(0,int(N/LL)):
    for j in range(0,LL):
        x_comb.append(i)
x_comb = np.array(x_comb)

for i in range(1,31):

    J_L = (tau**2)*np.ones((LL,LL))
    V_i = J_L + (LS_RSE)*np.identity(LL)
    V_inverse = np.linalg.inv(V_i)
    V_big = V_inverse
    for j in range(1,int(KL)):
        V_big = sc.linalg.block_diag(V_big,V_inverse)
    
    detV = KL*(np.log(np.absolute(np.linalg.det(V_i))))
    p_group = np.zeros((m_group,11))
    
    x_big = np.zeros((1,6))
    y_big = np.array([])
    

    for k in range(0,11):
        a = s[k]
        if a<0:
            backforward = 2*(np.ones((m_group,1)))
        elif a>0:
            backforward = np.ones((m_group,1))
        elif a ==0:
            backforward = np.zeros((m_group,1))
        
        shift = abs(a)*(np.ones((m_group,1)))
        groups = np.zeros((m_group,2))
        count = 0
        for l in range(0,m_group):
            groups[l,0] = int(count)
            groups[l,1] = int((l+1)*timepoints)
            count = int(count+timepoints)

        x_switch = np.zeros((m_group*timepoints,6))
        for l in range(0,m_group):
            index1 = int(groups[l,0])
            index2 = int(groups[l,1])
            if backforward[l] ==0:
                x_switch[index1:index2,:]=x_matrix[index1:index2,:]
            if backforward[l]==1:
                position1 = groups[l,0]+shift[l]
                position2 = groups[l,1]
                position1 = int(np.squeeze(position1))
                position2 = int(np.squeeze(position2))
                piece1 = x_matrix[position1:position2,:]
                position1 = groups[l,0]
                position2 = groups[l,0]+shift[l]
                position1 = int(np.squeeze(position1))
                position2 = int(np.squeeze(position2))
                piece2=x_matrix[position1:position2,:]
                x_switch[index1:index2,:]=np.concatenate((piece1, piece2))
            if backforward[l] == 2:
                position1 = groups[l,1]-shift[l]
                position2 = groups[l,1]
                position1 = int(np.squeeze(position1))
                position2 = int(np.squeeze(position2))
                piece1 = x_matrix[position1:position2,:]
                position1 = groups[l,0]
                position2 = groups[l,1]-shift[l]
                position1 = int(np.squeeze(position1))
                position2 = int(np.squeeze(position2))
                piece2=x_matrix[position1:position2,:]
                x_switch[index1:index2,:]=np.concatenate((piece1, piece2))

        RESP = DATA['value']
        resp = RESP.to_numpy()
        for l in range(0,m_group):
            index1 = int(groups[l,0])
            index2 = int(groups[l,1])
            M = resp[index1:index2]-np.matmul(x_switch[index1:index2,0:6], np.transpose(beta))
            p_group[l,k] = (-0.5*timepoints*np.log(2*np.pi) -0.5*detV -0.5*np.matmul( np.matmul( M , V_big ), np.transpose(M)))

        x_big=np.concatenate((x_big,x_switch))
        y_big= np.concatenate((y_big,resp))

    x_big=x_big[1:len(x_big[:,0]),:]
    #y_big=y_big[1:len(y_big[0])] come back maybe im not sure what this achieves

    for ii in range(0,m_group):
        p_group[ii,:]=np.exp(p_group[ii,:]-max(p_group[ii,:]))/sum(np.exp(p_group[ii,:]-max(p_group[ii,:])))
    
    print(p_group)
    print("")
    weights1 = np.zeros((N_total,11))

    for kk in range(0,11):
        for ii in range(0,m_group):
            index1 = int(groups[ii,0])
            index2 = int(groups[ii,1])
            weights1[index1:index2,kk]=p_group[ii,kk]
    weights2 = weights1[:,0]

    for kk in range(1,11):
        weights2=np.concatenate((weights2,weights1[:,kk]))
    
    weights2 = weights2*11

    data_big = np.column_stack((y_big,x_big))
    data_big = np.column_stack((data_big,weights2))
    data_big=np.column_stack((data_big,x_comb))
    data_big = pd.DataFrame(data=data_big, columns= ["resp","x1","x2","x3","x4","x5","x6","weights2","xcomb"])

    #data_big["groups"]=np.arange(0,len(data_big)) 
    #data_big = data_big.groupby('xcomb')['resp'].mean()
    data_big['groups'] = range(0, len(data_big) )
    grouped_data = data_big.groupby('xcomb')['resp'].apply(list)
    data_big_grouped = grouped_data.reset_index()
    data_big_grouped.columns = ['xcomb', 'resp']

    data_big["weights2"]= weights2 + 0.00000000000001
    weights22 = 1 / data_big['weights2']
    weights22[weights22 == np.inf] = np.nan  # Replace infinite values with NaN
    weights22.fillna(0.0000000000000001, inplace=True)

    
    #mm = sm.WLS.from_formula("resp ~ -1 + x1 + x2 + x3 + x4 + x5 + x6", data=data_big, weights=weights22)
    #result = mm.fit(cov_type='nonrobust')

    mm = sm.WLS(data_big['resp'], data_big[['x1', 'x2', 'x3', 'x4', 'x5', 'x6']],  weights=weights22)
    result = mm.fit() #gave working
    
   # mm = sm.GLS(data_big['resp'], data_big[['x1', 'x2', 'x3', 'x4', 'x5', 'x6']], weights=weights22)

    #grouped_data = data_big.groupby('xcomb')['resp'].apply(list)
    #correlation = sm.cov_struct.CompoundSymmetry().from_formula("1", groups=grouped_data)

    #mm = mm.fit(cov_type='cluster', cov_kwds={'groups': data_big['xcomb'], 'cov_struct': correlation})

    #result = mm
    beta3 = result.params.values.reshape(-1, 1)
    rho_star = result.cov_params().iloc[0, 0] 
    sig_star = result.scale
    tau = np.sqrt(rho_star * sig_star)
    sigma = np.sqrt(sig_star * (1 - rho_star))

    mm3 = result


    #mm3 = sm.MixedLM(endog=weights22*data_big['resp'], exog=data_big[["x1","x2","x3","x4","x5","x6"]], groups=data_big['xcomb'], exog_re=None)
    #mm3 = mm3.fit()


	#data_big <- groupedData(resp ~ 1 | xcomb, data_big)
	#data_big[,"weights2"]=data_big[,"weights2"]+0.0000000000000001
    #mm=gls(resp ~ -1+x1+x2+x3+x4+x5+x6,data=data_big,correlation=corCompSymm(value=0.3,form=~1|xcomb),weights=~1/weights2,method="ML")
    #beta = np.zeros((6,1))
########



tot_switch=11*m_group
groups_switch=np.zeros((tot_switch,4))
counter=0
groups_switchcol2 = list(np.arange(1,19))
groups_switch[:,1]=np.array(11*groups_switchcol2)
counter_2=0
counter_3=0
for i in range(0,tot_switch):
    groups_switch[i,2]=counter
    groups_switch[i,3]=(i+1)*timepoints
    counter=counter+timepoints
    groups_switch[i,0]=counter_2+1
    counter_3=counter_3+1
    if counter_3==18:
        counter_2=counter_2+1
        counter_3=0

x_back=np.zeros((1,6))
p_group = np.round(p_group)

kk2=6*np.ones((m_group,1))
for i in range(0,m_group):
    for j in range(0,11):
        if int(p_group[i,j])==1:
            kk2[i]=j

for j in range(0,m_group):
    for i in range(0,int(np.prod(m_group*11))):
        if (groups_switch[i,0]==kk2[j] and groups_switch[i,1]==j+1):
            index1 = int(groups_switch[i,2])
            index2 = int(groups_switch[i,3])
            x_back=np.concatenate((x_back,x_big[index1:index2,:]))
            #print(x_back)
x_back=x_back[1:len(x_back[:,0]),:]


mm = sm.OLS(resp, sm.add_constant(x_back[:,0:6])).fit()
mm1=mm
sigma_squared= (mm1.resid.std(ddof=(x_back[:,0:6]).shape[1]))**2

beta1 = mm1.fittedvalues
RSE = sum(np.square(resp-beta1))/3240


data_big = np.column_stack((resp,x_back,x_comb[0:N_total]))
data_big = pd.DataFrame(data=data_big, columns= ["resp","x1","x2","x3","x4","x5","x6","x_comb"])
data_big.groupby(['x_comb']).mean()

    #data_big["groups"]=np.arange(0,len(data_big)) 
    #data_big = data_big.groupby('xcomb')['resp'].mean()

#grouped_data = data_big.groupby('xcomb')['resp'].apply(list)
#data_big_grouped = grouped_data.reset_index()
#data_big_grouped.columns = ['xcomb', 'resp']

mm3 = sm.MixedLM.from_formula("resp ~ -1 + x1 + x2 + x3 + x4 + x5 + x6", groups="x_comb", data=data_big ,re_formula="1")
result3 = mm3.fit(method='powell', maxiter=1000)
beta = result3.params.values.reshape(-1, 1)

print(result3.summary())

varcorr_matrix = result3.cov_re
#tau = np.sqrt(varcorr_matrix.iloc[0, 1])
#sigma = np.sqrt(varcorr_matrix.iloc[1, 1])

tau =1
sigma=1
J_L = (tau**2)*np.ones((LL,LL))
V_i = J_L + (sigma**2)*np.identity(LL)
V_inverse = np.linalg.inv(V_i)
V_big = V_inverse
for j in range(0,int(KL)-1):
    V_big = sc.linalg.block_diag(V_big,V_inverse)
V_big2 = V_big
for j in range(0,m_group-1):
    V_big2 = sc.linalg.block_diag(V_big2,V_big)

SE = np.linalg.inv(np.matmul(np.transpose(x_back),np.matmul(V_big2,x_back)))

mm0f = mm0.fittedvalues
mm1f = mm1.fittedvalues
mm3f = result3.fittedvalues
####Plots
StationNo = 15
Start = StationNo*180
End = (StationNo+1)*180
plt.figure(figsize=(12, 8))
plt.plot(np.arange(1,181), resp[Start:End], 'k.', label='observed')
plt.plot(np.arange(1,181), mm0f[Start:End], 'b-', label='least squares')
plt.plot(np.arange(1,181), mm1f[Start:End], 'r-', label='timeshift')
plt.plot(np.arange(1,181), mm3f[Start:End],  'g-', label='timeshift+amplitude')
plt.legend()
plt.xlabel('Time')
plt.ylabel('Displacement')