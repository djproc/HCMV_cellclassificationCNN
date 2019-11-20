#writing functions to be called from different notebooks, making the code easier to read
import os
import glob
import pandas as pd

#TODO creat a nice function that saves a file summarizing the data outputs as graphs from cellprofiler
#def DataCheckPlotting():

def ImportData_NUC_CYTO(cp_output_dir):
    #use glob to pull the common filename from the csv files, call this the experiment name 
    #import glob
    #import os
    exp_name = os.path.basename(glob.glob(cp_output_dir + '/*_Image.csv',recursive=True)[0])[:-10]
    print(f'EXPERIMENT NAME: {exp_name}')

    #create a folder to export data from the CNN
    cnn_export_folder = f"_IMGexportforCNN_{exp_name}"
    os.makedirs(cnn_export_folder, exist_ok=True)

    #generate the other cellprofiler output filenames
    nuc_csv =  f"{exp_name}_NUC_DAPI.csv"
    cyto_csv = f"{exp_name}_Cytoplasm.csv"
    image_csv = f"{exp_name}_Image.csv"
    print("\nIMPORTED AND MERGED THE FOLLOWING FILES:", nuc_csv, cyto_csv, image_csv, sep="\n - ")

    #import these files as datafames using pandas
    #import pandas as pd
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
    print(f'\nDETECTED NUCLEI: {df_combined_data.shape[0]:,.0f}')
    return(df_combined_data, df_image_url,exp_name);

def GenerateIDsandCoords(df_combined_data, Metadata_for_ID):
    df_combined_data.reset_index(inplace=True)
    df_combined_data.rename(columns={'ObjectNumber_NUC':'NUC_ID'}, inplace=True)
    df_combined_data.dropna(subset=['NUC_ID'], inplace=True)
    
    #conversion to string for concatenation
    
    ### CODE CAN BE IMPROVED ###
    df_combined_data["NUC_ID"] = df_combined_data["NUC_ID"].astype(float).astype(int).astype(str)
    df_combined_data["ImageNumber"] = df_combined_data["ImageNumber"].astype(float).astype(int).astype(str)
    df_combined_data['Metadata_date_NUC'] = df_combined_data['Metadata_date_NUC'].astype(float).astype(int).astype(str)
    df_combined_data['Metadata_experiment_NUC'] = df_combined_data['Metadata_experiment_NUC'].astype(float).astype(int).astype(str)
    
    
    # TODO: if this is a string and not a number, then make it an float then int then string. 
    df_combined_data['Metadata_for_ID'] = df_combined_data[Metadata_for_ID].astype(float).astype(int).astype(str)
    
    df_combined_data["Unique_ID"] = df_combined_data['Metadata_experiment_NUC'].str.cat(
        df_combined_data[['Metadata_date_NUC', 'Metadata_biorep_NUC','Metadata_for_ID','ImageNumber','NUC_ID']]
        , sep="_")
    print(f'EXAMPLE IDS: {df_combined_data["Unique_ID"][1]}')
    return(df_combined_data);

def IDsandCoords_IMGexport(df_merge_OUTPUT,C1,C2,C3,exp_name):
    
    df_merge_OUTPUT['NUC_x0'] = df_merge_OUTPUT[('Location_Center_X_NUC')]*8 #correcting for downscaling
    df_merge_OUTPUT['NUC_y0'] = df_merge_OUTPUT[('Location_Center_Y_NUC')]*8 #correcting for downscaling
    
    df_merge_OUTPUT['CYTO_x0'] = df_merge_OUTPUT[('Location_CenterMassIntensity_X_gB_small')]*8 #correcting for downscaling
    df_merge_OUTPUT['CYTO_y0'] = df_merge_OUTPUT[('Location_CenterMassIntensity_Y_gB_small')]*8 #correcting for downscaling
    
    df_merge_OUTPUT['URL_C1'] = df_merge_OUTPUT[C1] #this will be red in imageJ
    df_merge_OUTPUT['URL_C2'] = df_merge_OUTPUT[C2] #this will be green in imageJ
    df_merge_OUTPUT['URL_C3'] = df_merge_OUTPUT[C3] #this will be blue in imageJ
    
    df_merge_EXPORT= df_merge_OUTPUT[["Unique_ID", 'NUC_x0','NUC_y0','URL_C1','URL_C2','URL_C3','AreaShape_Orientation','CYTO_x0','CYTO_y0']]
    
    df_merge_EXPORT.to_csv(f'{exp_name}_IDsandCoords.csv', index=False)
    print(f'Exported Nucleus IDs and Coordinates into a csv named: {exp_name}_IDsandCoords.csv')
    
    
    


