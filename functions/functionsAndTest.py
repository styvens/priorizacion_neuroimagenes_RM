import numpy as np
from sklearn.decomposition import PCA
import json
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import entropy
from scipy.stats import mannwhitneyu
from scipy import stats

def imagCentralTend(data,measure):
    if(measure=='mean'):
        cal_right = np.mean(data['right'],axis=0)
        cal_left = np.mean(data['left'],axis=0)
    elif(measure=='median'):
        cal_right = np.median(data['right'],axis=0)
        cal_left = np.median(data['left'],axis=0)
    
    return(cal_right,cal_left)

def imagCov(data):
    covR = np.cov(data['right'])
    covL = np.cov(data['left'])
    
    return(covR,covL)


def CompNum(s,threshold):
    sum_eig = np.sum(s)
    percentage_variance = np.divide(s, sum_eig)
    sum_var = 0
    num_var = 0
    for i in np.arange(percentage_variance.shape[0]):
        if sum_var >= threshold:
            num_var = i
            break;

        sum_var += percentage_variance[i]

    return(num_var)

def imagPCA(data,components):
    '''
    data: left and right lists
    componets: list with two values, firs corespond to right and econd to left
    '''
    pca_rigth = PCA(n_components=components[0], svd_solver='full').fit(data['right'])
    pca_left = PCA(n_components=components[1], svd_solver='full').fit(data['left'])
    return(pca_rigth,pca_left)

def hemisCOmpar(right,left,order):
    nR = np.linalg.norm(right,axis=1,ord=order)
    nL = np.linalg.norm(left,axis=1,ord=order)
    diff = np.abs(np.subtract(nR,nL))
    return(nR,nL,diff)


class NumpyEncoder(json.JSONEncoder):
    # function taken from https://interviewbubble.com/typeerror-object-of-type-float32-is-not-json-serializable/
    """ Special json encoder for numpy types """
    def default(self, obj):
        if isinstance(obj, (np.int_, np.intc, np.intp, np.int8,
            np.int16, np.int32, np.int64, np.uint8,
            np.uint16, np.uint32, np.uint64)):
            return int(obj)
        elif isinstance(obj, (np.float_, np.float16, np.float32, 
            np.float64)):
            return float(obj)
        elif isinstance(obj,(np.ndarray,)): #### This is the fix
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)
    
    
    
def thersPH(data_healty,data_pathologic,W,nboost):
    '''
    data: left and right values
    W: list witn two values, first is the weight of mean and second is the weight of proportion
    '''
    nH = len(data_healty)
    nP = len(data_pathologic)
    PH = []; PP = []
    MAXHO = []; CHECKP0 = []; MAXH = []; CHECKP = [];
    res_H=[];res_P = [];res_H0=[];res_P0 = []
    
    for b in range(nboost):   
    
        rH = pd.unique(np.random.choice(nH,nH, replace=True))
        rP = pd.unique(np.random.choice(nP,nP, replace=True))
        dat_healty = [data_healty[i] for i in rH] 
        dat_pathol = [data_pathologic[i] for i in rP] 
        
        for i in range(len(dat_healty)):
            difH = np.abs(np.subtract(dat_healty[i]['left'],dat_healty[i]['right']))
            pH = len(difH[difH>0])/len(difH)
            calH0 = np.mean(difH[difH>0]) + pH
            calH = np.mean(difH[difH>0])*W[0] + pH*W[1]
            res_H0.append(calH0)
            res_H.append(calH)
            PH.append(pH)

        for j in range(len(dat_pathol)):
            difP = np.abs(np.subtract(dat_pathol[j]['left'],dat_pathol[j]['right']))
            pP = len(difP[difP>0])/len(difP)
            calP0 = np.mean(difP[difP>0]) + pP
            calP = np.mean(difP[difP>0])*W[0] + pP*W[1]
            res_P0.append(calP0)
            res_P.append(calP)
            PP.append(pP)

        maxH0 = np.max(res_H0)
        filtered_check0 = [number for number in res_P0 if number > maxH0]
        checkP0 = len(filtered_check0)/len(res_P0)
        MAXHO.append(maxH0)
        CHECKP0.append(checkP0)

        maxH = np.max(res_H)
        filtered_check = [number for number in res_P if number > maxH]
        checkP = len(filtered_check)/len(res_P)
        MAXH.append(maxH)
        CHECKP.append(checkP)
    

    cehck = {'No_weight':CHECKP0,'weight':CHECKP}
    trheshol = {'No_weight':MAXHO,'weight':MAXH}
    
    return(trheshol,cehck,PH,PP,res_H,res_P,res_H0,res_P0)


def permut(n):
    '''
    Generation of uniform random numbers between 0 and n-1 with replace, where n is the size of a vector.
    '''
    r = pd.unique(np.random.choice(n,n, replace=True))
    nr = len(r)
    return(nr)

def ranselect(x,y,n):
    '''
    Select and exchange of values between two vectors, given the positions
    '''
    

    r = pd.unique(np.random.choice(n,n, replace=True)).tolist() 
    xchan = x[r]
    ychan = y[r]
    x[r] = ychan
    y[r] = xchan
   
    return(x,y)

def rCal(x,y,order,calType):
    '''
    Calculation of the norm of two vectors, each separately and the ratio between these two values.
    '''
    d0 = np.linalg.norm(x,ord=order)
    d1 = np.linalg.norm(y,ord=order)
    if(calType=='rate'):
        cal = np.divide(d0,d1)
    elif(calType=='difference'):
        cal = np.subtract(d0,d1)
    res = (d0,d1,cal)
    return(res)

def rPermut(n_permuta,x_l,x_r,order,n,calType):
    '''
    Execute the permut, ranselect and rCal functions n_permuta times to obtain n dictionaries with values 
    '''
    
    Nx = [];Ny=[];Res=[]
    Results = {}
    
    for  j in np.arange(n_permuta):
        #r = permut(n)
        X, Y = ranselect(x_l,x_r,n)
        cal=rCal(X,Y,order,calType)
        Nx.append(cal[0])
        Ny.append(cal[1])
        Res.append(cal[2])
        
    Results['norm_x'] = Nx
    Results['norm_y'] = Ny
    Results['R'] = Res
    return(Results)


def ranselect2(x,y):
    '''
    Select and exchange of values between two vectors, given the positions
    '''
    
    pos = np.where(x!=0)[0].tolist()
    xnew = x[pos]
    ynew = y[pos]
    nnew = len(pos)

    r = pd.unique(np.random.choice(nnew,nnew, replace=True)).tolist() 
    xchan = xnew[r]
    ychan = ynew[r]
    xnew[r] = ychan
    ynew[r] = xchan
    return(xnew,ynew)

def rCal2(x,y,order,calType):
    '''
    Calculation of the norm of two vectors, each separately and the ratio between these two values.
    '''
    d0 = np.linalg.norm(x,ord=order)
    d1 = np.linalg.norm(y,ord=order)
    if(calType=='rate'):
        cal = np.divide(d0,d1)
    elif(calType=='difference'):
        cal = np.subtract(d0,d1)
    res = (d0,d1,cal)
    return(res)

def rPermut2(n_permuta,x_l,x_r,order,n,calType):
    '''
    Execute the permut, ranselect and rCal functions n_permuta times to obtain n dictionaries with values 
    of norms and ratios
    '''
    Nx = [];Ny=[];Res=[]
    Results = {}
    
    for  j in np.arange(n_permuta):
        xnew, ynew = ranselect2(x_l,x_r)
        cal=rCal2(xnew,ynew,order,calType)
        Nx.append(cal[0])
        Ny.append(cal[1])
        Res.append(cal[2])
        
    Results['norm_x'] = Nx
    Results['norm_y'] = Ny
    Results['R'] = Res
    return(Results)

def distTest(n_figure,x,y,data,type_patient,norm,col):
    """
    n_figre<np.arange  """
    f=plt.figure(figsize=(16,14))
    dim = str(x)+str(y)
    for i in np.arange(n_figure):
        #ax=plt.subplot(x, y, i+1)
        ax=f.add_subplot(x,y,i+1)
        
        sns.distplot(data[i]['R'], hist = True, kde = True,ax=ax,
                         kde_kws = {'linewidth': 3},
                         label = data[i]['label'],color=col)

        
def compareTest(data):
    res_mean=[];q005=[];q025=[];q05=[];q25=[];q5=[];q75=[];q95=[];q975=[];q995=[];Range=[]
    U0=[];Un=[];std=[]
    for i in np.arange(len(data)):
        comp_mean=np.mean(data[i]['R'])
        Q = np.quantile(data[i]['R'],[0.005,0.025,0.05,0.25,0.5,0.75,0.095,0.975,0.995])
        q005.append(Q[0]);q025.append(Q[1]);q05.append(Q[2]);q25.append(Q[3]);q5.append(Q[5]);q75.append(Q[5])
        q95.append(Q[6]);q975.append(Q[7]);q995.append(Q[7])
        u0 = np.min(data[i]['R'])
        un = np.max(data[i]['R'])
        sd = np.std(data[i]['R'])
                    
        rang = un - u0
        Range.append(rang)
        U0.append(u0)
        Un.append(un)            
        res_mean.append(comp_mean)
        std.append(sd)
        
    res={'d_mean':res_mean,'d_q005':q005,'d_q025':q025,'d_q05':q05,'d_q25':q25,'d_q5':q5,'d_q75':q75,
         'd_q95':q95,'d_q975':q975,'d_q995':q995,'range':Range,'U0':U0,'Un':Un}
    
    print('Mean ratio: ' ,np.round(np.mean(res_mean),4),
         '\nMedian ratio: ' ,np.round(np.median(q5),4),
        '\nStandar deviation ratio: ' ,np.round(np.std(Range),4),
         '\nMean range_ratio: ' ,np.round(np.mean(Range),4),
         '\nInterquratile ratio: ' ,np.round(np.mean(q75)-np.mean(q25),4),
          '\nU0: ' ,np.round(np.mean(U0),4),
          '\nUn: ' ,np.round(np.mean(Un),4),
         '\nPercentile Uo: ' ,np.quantile(U0,[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]),
         '\nPercentile Un: ' ,np.quantile(Un,[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]),
         '\nPercentile std: ' ,np.quantile(std,[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]))
    return(res)


def checkHCP(dat,threshold0,thresholdn,thresholdsd,reff):
    res = []
    for i in range(len(dat)):
        data=np.asanyarray(dat[i]['R'])
        p0 = len(data[data<=threshold0])/len(data)
        pn = len(data[data>=thresholdn])/len(data)
        sd = np.std(data)
        
        if((p0>reff or pn>reff) and (sd>thresholdsd)):
            clas = 1
        elif((p0>reff or pn>reff) and (sd<=thresholdsd)):
            clas = 2
        else:
            clas = 0
        res.append(clas)
    return(res)


def detectLabelHigh(data,low,upp):
    true_class = []; false_class = []; all_labels = []
    for i in range(len(data)):
        lab = res_permut_pathologic[i]['label'][0]
        all_labels.append(lab)
        val = data[i]
        if(val>=low and val<=upp):
            false_class.append(lab)
        else:
            true_class.append(lab)
    
    print('Correctly clasified',len(true_class),
         '\nIncorrectly clasified',len(false_class),
         '\nTrue percent classified: ',np.round(len(true_class)/(len(true_class)+len(false_class))*100,3))
    return(true_class,false_class,all_labels)

def getDataFromDirectoryBianca(directories,nx,ny,nz,nt):
    images_healthy = []
    images_pathologic = []
    
    for f in directories:
        nam = f[0:11]
        l_nam = nam + '_bianca_output_l_flipped.nii.gz'
        X={}
        for f1 in directories:
            r_nam = nam + '_bianca_output_r.nii.gz'
        
            if(l_nam == f and r_nam == f1): # procesing only left and right hemispheres
                x = get_data(path+'/'+ nam + '*_bianca_output_l_flipped.nii.gz').reshape(nx*ny*nz*nt)
                y = get_data(path+'/'+ nam + '*_bianca_output_r.nii.gz').reshape(nx*ny*nz*nt)
                X['left'] = x;X['right']=y
                label=df_labels.loc[df_labels['id1']==nam,'note'].values
                X['label']=label

                if(label=='normal'):
                    images_healthy.append(X)
                else:
                    images_pathologic.append(X)

    return(images_healthy,images_pathologic)

def checkHCPBianca(dat,threshold0,thresholdn,thresholdsd,reff):
    res = []
    for i in range(len(dat)):
        data=np.asanyarray(dat[i]['R'])
        sd = np.std(data)
        
        if((sd>thresholdsd)):
            clas = 1
        elif((sd<=thresholdsd)):
            clas = 0
        res.append(clas)
    return(res)






def entropyData(x,y,nextend):
    '''
    Select and exchange of values between two vectors, given the positions
    input:
            x, y: array 3-dimension
            nextend: number of pixels around each dimension 
    '''
    dim1,dim2,dim3, _ = x.shape
    ref = np.argwhere(x)
    X = []; Y = []
    
    for k in range(ref.shape[0]):

        if(ref[k][0]<nextend):
            li1=0;ls1=nextend*2-1
        elif(ref[k][0]>dim1-nextend-1):
            li1=dim1-nextend-1;ls1=dim1-1
        else:
            li1=ref[k][0]-nextend;ls1=ref[k][0]+nextend
            
        if(ref[k][1]<nextend):
            li2=0;ls2=nextend*2-1
        elif(ref[k][1]>dim2-nextend-1):
            li2=dim2-nextend-1;ls2=dim2-1
        else:
            li2=ref[k][1]-nextend; ls2=ref[k][1]+nextend
        
        if(ref[k][2]<nextend):
            li3=0;ls3=nextend*2-1
        elif(ref[k][2]>dim3-nextend-1):
            li3=dim3-nextend-1;ls3=dim3-1
        else:
            li3=ref[k][2]-nextend; ls3=ref[k][2]+nextend
        
        
        datx = x[li1:ls1,  li2:ls2,  li3:ls3]
        dimx, dimy, dimz ,_ = datx.shape
        datx = datx.reshape(dimx*dimy*dimz)
        daty = y[li1:ls1,  li2:ls2,  li3:ls3].reshape(dimx*dimy*dimz)
        
        X.append(datx)
        Y.append(daty)
        
    return(X,Y)

def entropyRpermut(x,y,quantiles,nextend):
    '''
    x,y: right and lefth data
    quantiles: list with probabilitys
    '''
    ENTROPY_x = [];ENTROPY_y = []
    
    X, Y = entropyData(x,y,nextend)
   
    for k in range(len(X)):
        datx = X[k]
        daty = Y[k]
        nentro = len(datx)
        
        datnewx = datx; datnewy = daty
        # permutation
        r = pd.unique(np.random.choice(nentro,nentro, replace=True)).tolist()
        datxchan = datnewx[r]
        datychan = datnewy[r]
        datnewx[r] = datychan
        datnewy[r] = datxchan

        # Entropy
        Q_x = np.quantile(datx,quantiles)
        nq_x = int(np.round(100/(len(Q_x)),0))
        values_x = np.repeat(Q_x, nq_x, axis=0)
        _, counts_elements_x = np.unique(values_x, return_counts=True)
        cal_x = entropy(counts_elements_x)
        
        Q_y = np.quantile(daty,quantiles)
        nq_y = int(np.round(100/(len(Q_y)),0))
        values_y = np.repeat(Q_y, nq_y, axis=0)
        _, counts_elements_y = np.unique(values_y, return_counts=True)
        cal_y = entropy(counts_elements_y)
        #Entropy.append(cal)
        
        ENTROPY_x.append(cal_x)
        ENTROPY_y.append(cal_y)
        
    return(ENTROPY_x,ENTROPY_y)


def imageCrossPermutEntropy1(data,quantiles,nextend):
    '''
    Execute the  entropyData and entropyRpermut
     input:
            permut, ranselect and rCal inputs
    output:
            Results: dictionary with entropy and test values.
    '''

    raw_entropy_x = [];mean_entropy_x=[]; raw_entropy_y =[];mean_entropy_y=[];Label=[]
    mannwhitney_stast=[];mannwhitney_p=[];t_test_stast=[];t_test_p=[]
    Results = {}
    
    for  i in range(len(data)):
        x = data[i]['right']
        y = data[i]['left']
        label = data[i]['label']
        
        Entropy_x, Entropy_y = entropyRpermut(x,y,quantiles,nextend)
        meanentropy_x = np.mean(Entropy_x)
        raw_entropy_x.append(Entropy_x)
        mean_entropy_x.append(meanentropy_x)
        
        meanentropy_y = np.mean(Entropy_y)
        raw_entropy_y.append(Entropy_y)
        mean_entropy_y.append(meanentropy_y)
        Label.append(label)
        
        U_stat, U_p =mannwhitneyu(Entropy_x, Entropy_y)
        t_stat, t_p = stats.ttest_ind(Entropy_x, Entropy_y, equal_var = False)
        mannwhitney_stast.append(U_stat);mannwhitney_p.append(U_p)
        t_test_stast.append(t_stat);t_test_p.append(t_p)
    
    Results['raw_entropy_x'] = raw_entropy_x
    Results['mean_entropy_x'] = mean_entropy_x
    Results['raw_entropy_y'] = raw_entropy_y
    Results['mean_entropy_y'] = mean_entropy_y
    Results['mannwhitney_stast'] = mannwhitney_stast
    Results['mannwhitney_p'] = mannwhitney_p
    Results['t_test_stast'] = t_test_stast
    Results['t_test_p'] = t_test_p
    Results['label'] = Label
    
    return(Results)

def HPLE(data_H,data_P,p_ref,test):
    TP=[];TN=[];FP=[];FN=[];Acc=[];Precision=[];Recall=[];F1_score=[]
    res={}
    if(test=='mannWhitney'):
        xH = np.asanyarray(data_H['mannwhitney_p'])
        xP = np.asanyarray(data_P['mannwhitney_p'])
    elif(test=='t_test'):
        xH = np.asanyarray(data_H['t_test_p'])
        xP = np.asanyarray(data_P['t_test_p'])
    for i in range(len(p_ref)):
        reff = p_ref[i]
        tp = len(xP[xP<reff])
        tn = len(xH[xH>=reff])
        fp = len(xH[xH<reff])
        fn = len(xP[xP>=reff])
        acc = np.round((tp+tn)/(tp+tn+fp+fn),2)
        precision = np.round(tp/(tp+fp),2)
        recall = np.round(tp/(tp+fn),2)
        f1_score = np.round(2*(recall*precision)/(recall+precision),2)
        
        TP.append(tp);TN.append(tn);FP.append(fp);FN.append(fn)
        Acc.append(acc);Precision.append(precision);Recall.append(recall);F1_score.append(f1_score)
    
    res['TP']=TP;res['TN']=TN;res['FP']=FP;res['FN']=FN
    res['Acc']=Acc;res['Precision']=Precision;res['Recall']=Recall;res['F1_score']=F1_score
    return(res)


def imagCentralTend_eigen(data,measure):
    if(measure=='mean'):
        cal_right = np.mean(data['right'],axis=0)
        cal_left = np.mean(data['left'],axis=0)
    elif(measure=='median'):
        cal_right = np.median(data['right'],axis=0)
        cal_left = np.median(data['left'],axis=0)
    
    return(cal_right,cal_left)

def imagCov_eigen(data,space):
    if(space=='hemispheres'):
        covR = np.cov(data['right'])
        covL = np.cov(data['left'])
        return(covR,covL)
    elif(space=='all'):
        cov = np.cov(data['data'])
        return(cov)

def imagPCA_eigen(data,components,space):
    '''
    data: left and right lists
    componets: list with two values, first corespond to right and second to left
    '''
    if(space=='hemispheres'):       
        pca_rigth = PCA(n_components=components[0], svd_solver='full').fit(data['right'])
        pca_left = PCA(n_components=components[1], svd_solver='full').fit(data['left'])
        return(pca_rigth,pca_left)
    elif(space=='all'):
        pca = PCA(n_components=components[0], svd_solver='full').fit(data['data'])
        return(pca)
    
