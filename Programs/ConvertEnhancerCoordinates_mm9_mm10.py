import pandas as pd
import pyliftover

def Converter(genome_build_in='mm9', genome_build_out='mm10'):
    return pyliftover.LiftOver(genome_build_in, genome_build_out)

def convertEnhancerDF(data, Converter, feature_name=None):
    cols=list(data.columns)
    if 'name' in data.columns:
        cols.remove('name')
    if len(cols) > 6:
        cols=cols[:6]
    cvrt_data = pd.DataFrame(columns=cols)
    cvrt_data.index.name='name'
    cvrt_idx = 1
    for idx in range(len(data)):
        cvrt_start = Converter.convert_coordinate(data.loc[idx]['chr'], data.loc[idx]['start'])
        cvrt_end = Converter.convert_coordinate(data.loc[idx]['chr'], data.loc[idx]['end'])
        if (cvrt_start is not None) and (cvrt_end is not None):
            if len(cvrt_start)==1 and len(cvrt_end)==1:                  # only one converted genomic region should exist
                if cvrt_start[0][0] == cvrt_end[0][0]:                  # feature must have both ends on the same chromosome
                    chrom = cvrt_start[0][0]
                    if cvrt_start[0][1] > cvrt_end[0][1]:
                        start = cvrt_end[0][1]
                        end = cvrt_start[0][1]
                    else:
                        start = cvrt_start[0][1]
                        end = cvrt_end[0][1]
                    if end > start:                                     # only take features that have a length greater than 0
                        if feature_name is not None:
                            index_name = feature_name + str(cvrt_idx)
                        else:
                            if 'name' in data.columns:
                                index_name=data.loc[idx]['name']
                            else:
                                index_name = cvrt_idx
                        row = [chrom, start, end]
                        if len(cols) > 3:
                            for col in cols[3:]:                        # transfer some information from columns
                                row.append(data.loc[idx][col])          # append new row with enhancer information, index is enhancer name
                        cvrt_data.loc[index_name]=row
                        cvrt_idx += 1
    return cvrt_data

DFmm9 = pd.DataFrame.from_csv("/media/linux/TOSHIBA EXT/CHIPSeq_Asun2019/CSV_feature_files/GSE117034_HIRA_Flag_vs_Input_peaks.csv", index_col=False)
DFmm10 = convertEnhancerDF(DFmm9, Converter('mm9', 'mm10'), feature_name='Hira_')
DFmm10.to_csv("/home/linux/Asun_ChIPseq/genomic_coordinates_features/ES_Hira_mm10.csv")

DFmm9 = pd.DataFrame.from_csv("/media/linux/TOSHIBA EXT/CHIPSeq_Asun2019/CSV_feature_files/GSE117034_H33_HA_vs_Input_peaks.csv", index_col=False)
DFmm10 = convertEnhancerDF(DFmm9, Converter('mm9', 'mm10'), feature_name='H3.3_')
DFmm10.to_csv("/home/linux/Asun_ChIPseq/genomic_coordinates_features/ES_H33_HA_mm10.csv")

DFmm9 = pd.DataFrame.from_csv("/media/linux/TOSHIBA EXT/CHIPSeq_Asun2019/CSV_feature_files/GSE117034_UBN2_Flag_vs_Input_peaks.csv", index_col=False)
DFmm10 = convertEnhancerDF(DFmm9, Converter('mm9', 'mm10'), feature_name='Ubn2_')
DFmm10.to_csv("/home/linux/Asun_ChIPseq/genomic_coordinates_features/ES_Ubn2_mm10.csv")

DFmm9 = pd.DataFrame.from_csv("/media/linux/TOSHIBA EXT/CHIPSeq_Asun2019/CSV_feature_files/GSE117034_UBN1_Flag_vs_Input_peaks.csv", index_col=False)
DFmm10 = convertEnhancerDF(DFmm9, Converter('mm9', 'mm10'), feature_name='Ubn1_')
DFmm10.to_csv("/home/linux/Asun_ChIPseq/genomic_coordinates_features/ES_Ubn1_mm10.csv")


"""
enh_mm9 = pd.DataFrame.from_csv("/home/linux/Asun_ChIPseq/genomic_coordinates_features/ES_Enhancers_mm9.csv", index_col=False)
enh_mm10 = convertEnhancerDF(enh_mm9, Converter('mm9', 'mm10'), feature_name='enh')
enh_mm10.to_csv("/home/linux/Asun_ChIPseq/genomic_coordinates_features/ES_Enhancers_mm10.csv")

prL_mm9=pd.DataFrame.from_csv("/home/linux/Asun_ChIPseq/genomic_coordinates_features/ES_PromoterLike_mm9.csv", index_col=False) 
prL_mm10 = convertEnhancerDF(prL_mm9, Converter('mm9', 'mm10'), feature_name='prL')
prL_mm10.to_csv("/home/linux/Asun_ChIPseq/genomic_coordinates_features/ES_PromoterLike_mm10.csv")

Xacc_mm9=pd.DataFrame.from_csv("/home/linux/Asun_ChIPseq/genomic_coordinates_features/TS_Xchr_accessibility_mm9.csv", index_col=False) 
Xacc_mm10 = convertEnhancerDF(Xacc_mm9, Converter('mm9', 'mm10'))
Xacc_mm10.to_csv("/home/linux/Asun_ChIPseq/genomic_coordinates_features/TS_Xchr_accessibility_mm10.csv")


Ezh2_d0_mm9 = pd.DataFrame.from_csv("/home/linux/Asun_ChIPseq/genomic_coordinates_features/JT_Lee_chrX_Ezh2_d0_strong_mm9.csv", index_col=False)
Ezh2_d0_mm10 = convertEnhancerDF(Ezh2_d0_mm9, Converter('mm9', 'mm10'))
Ezh2_d0_mm10.to_csv("/home/linux/Asun_ChIPseq/genomic_coordinates_features/Ezh2_d0_mm10.csv")

Ezh2_d7_mm9 = pd.DataFrame.from_csv("/home/linux/Asun_ChIPseq/genomic_coordinates_features/JT_Lee_chrX_Ezh2_d7_strong_mm9.csv", index_col=False)
Ezh2_d7_mm10 = convertEnhancerDF(Ezh2_d7_mm9, Converter('mm9', 'mm10'))
Ezh2_d7_mm10.to_csv("/home/linux/Asun_ChIPseq/genomic_coordinates_features/Ezh2_d7_mm10.csv")

Ezh2_d0_moderate_mm9 = pd.DataFrame.from_csv("/home/linux/Asun_ChIPseq/genomic_coordinates_features/JT_Lee_chrX_Ezh2_d0_moderate_mm9.csv", index_col=False)
Ezh2_d0_moderate_mm10 = convertEnhancerDF(Ezh2_d0_moderate_mm9, Converter('mm9', 'mm10'))
Ezh2_d0_moderate_mm10.to_csv("/home/linux/Asun_ChIPseq/genomic_coordinates_features/Ezh2_d0_moderate_mm10.csv")
"""
