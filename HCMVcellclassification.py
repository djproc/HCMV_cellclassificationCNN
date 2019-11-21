#writing functions to be called from different notebooks, making the code easier to read
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import shutil
%matplotlib inline
import matplotlib.pyplot as plt
import pandas as pd
import cv2
from tqdm import tqdm


def DataCheck(df_combined_data, group_data_by, display_x_axis, display_y_axis):
    fig, ax =plt.subplots(1,3, figsize=(20,4))

    sns.violinplot(x=group_data_by, y='AreaShape_Area', data=df_combined_data.reset_index(), size=1,ax=ax[0])
    sns.violinplot(x=group_data_by, y=display_y_axis, data=df_combined_data.reset_index(), size=1,ax=ax[1])
    sns.scatterplot(x=display_x_axis, y=display_y_axis, data=df_combined_data.reset_index(), hue=group_data_by, size=0.0001, alpha=0.01,ax=ax[2],legend=False)
    fig.show()
    #TODO: RENDER A FEW CELLS HERE ALSO, RANDOMLY SO WE GET AN IDEA OF DATA QUALITY BEFORE EXPORT

def ImportData_NUC_CYTO():
    cp_output_dir = "_CellProfiler"
    os.makedirs(cp_output_dir, exist_ok=True)
    #print(f'CREATED FOLDER NAMED: {cp_output_dir} \nFOLDER LOCATED AT: {os.getcwd()}')
    
    #TODO: get .csv files and add to this folder, also add any .cppipe file and any file with the same name
    exp_name = os.path.basename(glob.glob('*_Image.csv',recursive=True)[0])[:-10] #_Image.csv files will always be saved
    print(f'EXPERIMENT NAME: {exp_name}')
    
    #move files with exp_name in them to the cp_output_dir
    exp_files = glob.glob(f'{exp_name}*',recursive=True)
    for f in exp_files:
        shutil.move(f, os.path.join(cp_output_dir,f))
        
    #generate the other cellprofiler output filenames
    nuc_csv =  f"{exp_name}_NUC_DAPI.csv"
    cyto_csv = f"{exp_name}_Cytoplasm.csv"
    image_csv = f"{exp_name}_Image.csv"
    print("\nIMPORTED AND MERGED THE FOLLOWING FILES:", nuc_csv, cyto_csv, image_csv, sep="\n - ")
    
    #import these files as datafames using pandas
    #nucleus data
    df_nuc = pd.read_csv(os.path.join(cp_output_dir, nuc_csv), na_filter=True)
    df_nuc.set_index("ImageNumber", inplace=True)
    #cytoplasm data
    df_cyto = pd.read_csv(os.path.join(cp_output_dir, cyto_csv), na_filter=True)
    df_cyto.set_index("ImageNumber", inplace=True)
    #image info
    df_image = pd.read_csv(os.path.join(cp_output_dir, image_csv), na_filter=True)
    df_image.set_index("ImageNumber", inplace=True)
    #then extract only the image urls from this 
    df_image_url = df_image.filter(regex=r'^URL_', axis=1) #this will select any columns starting with "URL_"
    
    #combine these dataframes together
    #merge nucleus data with urls
    df_combined_data = df_nuc.merge(df_image_url, left_on='ImageNumber', right_on='ImageNumber', how='outer')
    #merge this with cytoplasm data and differentiate columns from the two datasets as "_NUC" and "_CYTO"
    df_combined_data = df_combined_data.merge(df_cyto, left_on=["ImageNumber", "ObjectNumber"], right_on=["ImageNumber", "Parent_NUC_DAPI"], how="outer", suffixes=('_NUC', '_CYTO'))

    #we can also just look at the raw number of rows in the dataframe to see how many nuclei we've identified
    df_combined_data.reset_index(inplace=True)
    df_combined_data.rename(columns={'ObjectNumber_NUC':'NUC_ID'}, inplace=True)
    df_combined_data.dropna(subset=['NUC_ID'], inplace=True)
    print(f'\nDETECTED NUCLEI: {df_combined_data.shape[0]:,.0f}')
    return(df_combined_data, df_image_url, exp_name);


def GenerateIDs_IMGexport(df_combined_data, group_data_by,C1,C2,C3,exp_name, img_size):
    
    #conversion to string for concatenation
    ### CODE CAN BE IMPROVED ###
    df_combined_data["NUC_ID"] = df_combined_data["NUC_ID"].astype(float).astype(int).astype(str)
    df_combined_data["ImageNumber"] = df_combined_data["ImageNumber"].astype(float).astype(int).astype(str)
    df_combined_data['Metadata_date_NUC'] = df_combined_data['Metadata_date_NUC'].astype(float).astype(int).astype(str)
    df_combined_data['Metadata_experiment_NUC'] = df_combined_data['Metadata_experiment_NUC'].astype(float).astype(int).astype(str) 
    
    # TODO: if this is a string and not a number, then make it an float then int then string. 
    df_combined_data['Metadata_for_ID'] = df_combined_data[group_data_by].astype(str)
    
    # If there is no biorep value, set it to 'repX'
    if 'Metadata_biorep_NUC' not in df_combined_data:
        df_combined_data['Metadata_biorep_NUC'] = 'repX' 
    
    df_combined_data["Unique_ID"] = df_combined_data['Metadata_experiment_NUC'].str.cat(df_combined_data[['Metadata_date_NUC', 
                                                                                                          'Metadata_biorep_NUC',
                                                                                                          'Metadata_for_ID',
                                                                                                          'ImageNumber','NUC_ID']
                                                                                                        ], sep="_")
    print(f'EXAMPLE IDS: {df_combined_data["Unique_ID"][1]}')
    
    df_combined_data['NUC_x0'] = df_combined_data[('Location_Center_X_NUC')]*8 #correcting for downscaling
    df_combined_data['NUC_y0'] = df_combined_data[('Location_Center_Y_NUC')]*8 #correcting for downscaling
    
    df_combined_data['CYTO_x0'] = df_combined_data[('Location_CenterMassIntensity_X_gB_small')]*8 #correcting for downscaling
    df_combined_data['CYTO_y0'] = df_combined_data[('Location_CenterMassIntensity_Y_gB_small')]*8 #correcting for downscaling
    
    df_combined_data['URL_C1'] = df_combined_data[C1] #this will be red in imageJ
    df_combined_data['URL_C2'] = df_combined_data[C2] #this will be green in imageJ
    df_combined_data['URL_C3'] = df_combined_data[C3] #this will be blue in imageJ
    
    df_display = df_combined_data[["Unique_ID",
                                   'NUC_x0',
                                   'NUC_y0',
                                   'URL_C1',
                                   'URL_C2',
                                   'URL_C3',
                                   'AreaShape_Orientation',
                                   'CYTO_x0',
                                   'CYTO_y0']]
        
    #create a folder to export data from the CNN
    cnn_export_folder = f"_IMGexportforCNN_{exp_name}"
    os.makedirs(cnn_export_folder, exist_ok=True)

    #just in case this crashes half way through, if you re-run this function this checks what you have already exported and continues from there
    files = [os.path.splitext(filename)[0] for filename in os.listdir(cnn_export_folder)]
    for index, row in tqdm(df_display.iterrows(), total=df_display.shape[0]):
        if df_display['Unique_ID'][index] not in files:
    
            img_border = int(img_size/2)
            
            ####TODO: this could all be far less memory intensive###
            red = cv2.imread(df_display['URL_C1'][index],-1)
            green = cv2.imread(df_display['URL_C2'][index],-1)
            blue = cv2.imread(df_display['URL_C3'][index],-1)

            red_8bit = (red/32).astype('uint8')
            green_8bit = (green/32).astype('uint8')
            blue_8bit = (blue/32).astype('uint8') * 0 # this makes it black 
            
            red = None
            green = None
            blue = None

            merged_channels = cv2.merge((red_8bit,green_8bit,blue_8bit))

            red_8bit = None
            green_8bit = None
            blue_8bit = None

            merged_channels = cv2.copyMakeBorder(merged_channels, img_border, img_border, img_border, img_border, cv2.BORDER_CONSTANT, value=0) 

            x0 = df_display['NUC_x0'][index]
            y0 = df_display['NUC_y0'][index]

            fig = plt.figure(frameon=False)
            fig.set_size_inches(2,2)
            ax = plt.Axes(fig, [0., 0., 1., 1.])
            ax.set_axis_off()
            fig.add_axes(ax)
            plt.xlim(x0,x0+img_size)
            plt.ylim(y0,y0+img_size)
            ax.imshow(merged_channels)

            filename = os.path.join(f'{cnn_export_folder}',(df_display['Unique_ID'][index] + '.jpg'))
            fig.savefig(filename, dpi=300, bbox_inches='tight', pad_inches=0)

            merged_channels = None

            fig.clf()
            plt.close(fig)
            plt.close('all')
         

        


def URLFIX(df, find_string, replace_string, columns):
    """Use this to fix path URL issues between operating systems"""
    for col in columns:
        df[col] = df[col].str.replace(find_string,replace_string)
        print(df[col][1])
    return(df)

def colSWAP(df, col_name_A, col_name_B):
    """If  you've named your metadata incrorrectly in CellProfiler, you can easily swap dataframe column names with this function"""
    df.rename(columns={col_name_A:col_name_B, col_name_B:col_name_A})
    return(df)

# TODO: linescan to analysis functions
"""
def LoadIDsfromCNN(exp_name, date_var, time_var, folder_name_1, folder_name_2):
    #exp_var = "137_20190805_SUN2quant_rerpCGH_v1"
    #date_var = 20190817
    #time_var = "1145"
    #### make this a drop down box ####
    #MOCK_confidece_list = [99,90]
    #TB_confidence_list = [99,90]

    FOLDER_1_files = []
    for conf in FOLDER_1_confidece_list:
        FOLDER_1_files_single_folder = os.listdir(f'{PATH}\\CNNpredictions_{date_var}_{time_var}\\{conf}confidence\\{folder_name_1}')
        print(f'Loading {len(MOCK_files_single_folder)} files for MOCK {conf}confidence')
        MOCK_files = MOCK_files + MOCK_files_single_folder
    print(f'Loaded {len(MOCK_files)} files total for TB')

    TB_files = []
    for conf in TB_confidence_list:
        TB_files_single_folder = os.listdir(f'{PATH}\\CNNpredictions_{date_var}_{time_var}\\{conf}confidence\\TB_perfect')
        print(f'Loading {len(TB_files_single_folder)} files for TB {conf}confidence')
        TB_files = TB_files + TB_files_single_folder
    print(f'Loaded {len(TB_files)} files total for TB')

    #convert these lists to a dataframe
    df_TB_ID = pd.DataFrame({'Unique_ID':TB_files})
    df_MOCK_ID = pd.DataFrame({'Unique_ID':MOCK_files})

    #strip "_RGB.jpg" from filename to match UNIQUE ID
    df_TB_ID['Unique_ID'].replace(regex=True,inplace=True,to_replace='_RGB.jpg',value='')
    df_MOCK_ID['Unique_ID'].replace(regex=True,inplace=True,to_replace='_RGB.jpg',value='')

    #filter mock for only in mock sample (there are obviously some uninfected cells in the TB infected samples)
    df_MOCK_ID = df_MOCK_ID[df_MOCK_ID['Unique_ID'].str.contains('MOCK')]
    print(df_MOCK_ID.shape)
    #filter TB for only in TB sample (just in case)
    df_TB_ID = df_TB_ID[df_TB_ID['Unique_ID'].str.contains('_96hpi_')]
    print(df_TB_ID.shape)

    # IMPORT DATAFRAME WHERE WE CREATED THE IDs and COORDS FOR THIS DATASET
    df_coords = pd.read_csv(f'{exp_var}_IDsandCoords.csv')
    #df_FULL = pd.read_csv("117_20181023_SUN2_v4_IDsandCoords_FULL.csv") #create a FULL version of the dataframe in the Generate IDs and coords notebook

    # Merge coords with lists from CNN sorting for TB
    df_TB_CNNsorted99_coords = df_coords.merge(df_TB_ID, left_on='Unique_ID', right_on='Unique_ID', how='right')
    df_MOCK_CNNsorted99_coords = df_coords.merge(df_MOCK_ID, left_on='Unique_ID', right_on='Unique_ID', how='right')
"""

# TODO: merging and graphing functions


