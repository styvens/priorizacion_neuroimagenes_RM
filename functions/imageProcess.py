import numpy as np
from nilearn import image
from nilearn.image import get_data


def getDataFromDirectory(directories,space,nx,ny,nz,centralice):
    images={}
    left=[]; right=[];labels=[]
    if(space=='mask'):
        for directory in directories:
            x = get_data(directory + '*_mask_l_flipped.nii.gz').reshape(nx*ny*nz)
            y = get_data(directory + '*_mask_r.nii.gz').reshape(nx*ny*nz)
            left.append(x);right.append(y)
            start = directory.find('/dicom/')
            end = len(directory)-1
            file=directory[start+7:end]
            label=df_labels.loc[df_labels['id1']==file,'note']
            labels.append(label.values)
            
    elif(space=='mni'):
        for directory in directories:
            x = get_data(directory + '*_mni-2mm_l_flipped.nii.gz').reshape(nx*ny*nz)
            y = get_data(directory + '*_mni-2mm_r.nii.gz').reshape(nx*ny*nz)
            left.append(x);right.append(y)
            start = directory.find('/dicom/')
            end = len(directory)-1
            file=directory[start+7:end]
            label=df_labels.loc[df_labels['id1']==file,'note']
            labels.append(label.values)
    
    # == Centralice image data
    if(centralice=='mean'):
        left = left-np.mean(left,axis=0)
        right = left-np.mean(right,axis=0)
    elif(centralice=='median'):
        left = left-np.median(left,axis=0)
        right = left-np.median(right,axis=0)
    elif(centralice==None):
        left = left
        rigth = right
    
    images['left'] = left;images['right'] = right; images['label'] = labels
    
    return images


def getDataFromDirectoryP(directories,space,nx,ny,nz):
    images=[]
    if(space=='mask'):
        for directory in directories:
            X={}
            x = get_data(directory + '*_mask_l_flipped.nii.gz').reshape(nx*ny*nz)
            y = get_data(directory + '*_mask_r.nii.gz').reshape(nx*ny*nz)
            X['left'] = x;X['right']=y
            start = directory.find('/dicom/')
            end = len(directory)-1
            file=directory[start+7:end]
            label=df_labels.loc[df_labels['id1']==file,'note']
            X['label']=label.values
            images.append(X)
    elif(space=='mni'):
        for directory in directories:
            X={}
            x = get_data(directory + '*_mni-2mm_l_flipped.nii.gz').reshape(nx*ny*nz)
            y = get_data(directory + '*_mni-2mm_r.nii.gz').reshape(nx*ny*nz)
            X['left'] = x;X['right']=y
            start = directory.find('/dicom/')
            end = len(directory)-1
            file=directory[start+7:end]
            label=df_labels.loc[df_labels['id1']==file,'note']
            X['label']=label.values
            images.append(X)
    return images

def getDataFromDirectoryHPLE(directories,space,nx,ny,nz):
    images=[]
    if(space=='mask'):
        for directory in directories:
            X={}
            x = get_data(directory + '*_mask_l_flipped.nii.gz')
            y = get_data(directory + '*_mask_r.nii.gz')
            X['left'] = x;X['right']=y
            start = directory.find('/dicom/')
            end = len(directory)-1
            file=directory[start+7:end]
            label=df_labels.loc[df_labels['id1']==file,'note']
            X['label']=label.values
            images.append(X)
    elif(space=='mni'):
        for directory in directories:
            X={}
            x = get_data(directory + '*_mni-2mm_l_flipped.nii.gz').reshape(nx,ny,nz)
            y = get_data(directory + '*_mni-2mm_r.nii.gz').reshape(nx,ny,nz)
            X['left'] = x;X['right']=y
            start = directory.find('/dicom/')
            end = len(directory)-1
            file=directory[start+7:end]
            label=df_labels.loc[df_labels['id1']==file,'note']
            X['label']=label.values
            images.append(X)
    return images

def getDataFromDirectory_NCCB_eigen(df_labels,directories,space,nx,ny,nz,centralice):
    images={}
    left=[]; right=[];labels=[];data=[]
    if(space=='hemispheres'):
        for directory in directories:
            x = get_data(directory + '_bianca_output_NCC_fl.nii.gz').reshape(nx,ny,nz)
            y = get_data(directory + '_bianca_output_NCC_r.nii.gz').reshape(nx,ny,nz)
            left.append(x);right.append(y)
            label=df_labels.loc[df_labels['path']==directory,'note'].values
            labels.append(label)
            
    elif(space=='all'):
        for directory in directories:
            x = get_data(directory + '_bianca_output_NCC_map_nt.nii.gz').reshape(nx,ny,nz)
            data.append(x)
            label=df_labels.loc[df_labels['path']==directory,'note'].values
            labels.append(label)
    
    # == Centralice image data
    if(centralice=='mean' and space=='hemispheres'):
        m_l=np.mean(left,axis=0); m_r=np.mean(right,axis=0)
        left = np.subtract(left,m_l); right = np.subtract(right,m_r)
    elif(centralice=='median' and space=='hemispheres'):
        me_l=np.median(left,axis=0); me_r=np.median(right,axis=0)
        left = np.subtract(left,me_l); right = np.subtract(right,me_r)
    elif(centralice=='mean' and space=='all'):
        m = np.mean(data,axis=0)
        data = np.subtract(data,m)
    elif(centralice=='median' and space=='all'):
        me=np.median(data,axis=0)
        data = np.subtract(data,me)
        
    if(space=='hemispheres'):
        left_new = [];right_new = []
        for l in range(len(left)):
            xnew = left[l].reshape(nx*ny*nz)
            left_new.append(xnew)
        
        for r in range(len(right)):
            ynew = right[r].reshape(nx*ny*nz)
            right_new.append(ynew)

        images['left'] = left_new; images['right'] = right_new; images['label'] = labels
    
    elif(space=='all'):
        data_new = []
        for i in range(len(data)):
            dnew = data[i].reshape(nx*ny*nz)
            data_new.append(dnew)
        
        images['data'] = data_new; images['label'] = labels
    
    return images

def getDataFromDirectory_NCC_eigen(directories,space,nx,ny,nz,centralice):
    images={}
    left=[]; right=[];labels=[];data=[]
    if(space=='hemispheres'):
        for directory in directories:
            x = get_data(directory + 'NCC_fl.nii.gz').reshape(nx,ny,nz)
            y = get_data(directory + 'NCC_r.nii.gz').reshape(nx,ny,nz)
            left.append(x);right.append(y)
            start = directory.find('/dicom/')
            end = len(directory)-1
            file=directory[start+7:end]
            label=df_labels.loc[df_labels['id1']==file,'note']
            labels.append(label.values)
            
    elif(space=='all'):
        for directory in directories:
            x = get_data(directory + 'NCC_map_nt.nii.gz').reshape(nx,ny,nz)
            data.append(x)
            start = directory.find('/dicom/')
            end = len(directory)-1
            file=directory[start+7:end]
            label=df_labels.loc[df_labels['id1']==file,'note']
            labels.append(label.values)
    
    # == Centralice image data
    if(centralice=='mean' and space=='hemispheres'):
        m_l=np.mean(left,axis=0); m_r=np.mean(right,axis=0)
        left = np.subtract(left,m_l); right = np.subtract(right,m_r)
        images['center_image']=m_l
    elif(centralice=='median' and space=='hemispheres'):
        me_l=np.median(left,axis=0); me_r=np.median(right,axis=0)
        left = np.subtract(left,me_l); right = np.subtract(right,me_r)
        images['center_image']=me_l
    elif(centralice=='mean' and space=='all'):
        m = np.mean(data,axis=0)
        data = np.subtract(data,m)
        images['center_image']=m
    elif(centralice=='median' and space=='all'):
        me=np.median(data,axis=0)
        data = np.subtract(data,me)
        images['center_image']=me
        
    if(space=='hemispheres'):
        left_new = [];right_new = []
        for l in range(len(left)):
            xnew = left[l].reshape(nx*ny*nz)
            left_new.append(xnew)
        
        for r in range(len(right)):
            ynew = right[r].reshape(nx*ny*nz)
            right_new.append(ynew)

        images['left'] = left_new; images['right'] = right_new; images['label'] = labels
    
    elif(space=='all'):
        data_new = []
        for i in range(len(data)):
            dnew = data[i].reshape(nx*ny*nz)
            data_new.append(dnew)
        
        images['data'] = data_new; images['label'] = labels
    
    return images
