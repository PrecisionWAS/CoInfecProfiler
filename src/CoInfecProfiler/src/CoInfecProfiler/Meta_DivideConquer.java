package CoInfecProfiler;

import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.time.format.DateTimeFormatter;
import java.util.Properties;

public class Meta_DivideConquer {
	public int reconstruction_start;
	public int reconstruction_end;
	public int region_length;
	String out_dir;
	public int [] region_start;
	public int [] region_end;
	public int [] tiling_region_start;
	public int [] tiling_region_end;
	String project_name;
	public String reference;
	public String SNV_cutoff;
	public String min_mapping_qual;
	public String min_read_length;
	public String max_insert_length;
	public String sequence_err;
	public String MEC_improvement_cutoff;
	public int num_threads;
	public String aBayesQR;
	public String [] cmds;
	public int thread_index;
	public int Max_L0L1_Regional_Haps;
	public double min_hap_freq;
	public double mismatch_tolerance;
	public String tool; 
	public double sum_weights;
	public double maf_weights;
	public double ld_weights;
	public String R;
	public int nGamma;
	int maxSuppSize;
	public double regression_gamma_min;
    public double regression_gamma_max;
    public double regression_lambda;
	
	
	final DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy/MM/dd HH:mm:ss");
	
	public Meta_DivideConquer(String parameter_file, int len) throws IOException, InterruptedException {
		
		InputStream is = new FileInputStream(parameter_file);
        Properties prop = new Properties();
        prop.load(is);
        this.out_dir= prop.getProperty("Output_Path")+"/" ;
        this.reconstruction_start = 1;
        this.reconstruction_end = len;
        this.region_length = len+1;
        this.project_name= prop.getProperty("Proj_Name");
        this.reference = prop.getProperty("Reference_Seq") ;
        this.SNV_cutoff = prop.getProperty("SNV_Cutoff") ;
        this.min_mapping_qual = prop.getProperty("Min_Mapping_Qual") ;
        this.min_read_length = prop.getProperty("Min_Read_Length") ;

        this.sum_weights = Double.parseDouble(prop.getProperty("Regression_One_Vector_Weight"));
        this.maf_weights = Double.parseDouble(prop.getProperty("Regression_Hap_MAF_Weight"));
        this.ld_weights = Double.parseDouble(prop.getProperty("Regression_Hap_LD_Weight"));
        this.R= prop.getProperty("Rscript_path") ;
        this.maxSuppSize= Integer.parseInt(prop.getProperty("Maximum_Haps_R"));
        this.nGamma= Integer.parseInt(
                prop.getProperty("Regression_n_Gamma"));
        
        this.regression_gamma_max=Double.parseDouble
        		(prop.getProperty("Regression_Gamma_Max"));
        this.regression_gamma_min=Double.parseDouble
        		(prop.getProperty("Regression_Gamma_Min"));
        
        this.regression_lambda = Double.parseDouble
        		(prop.getProperty("Regression_Lambda"));
        is.close();
        int num_regions = (this.reconstruction_end-this.reconstruction_start+1) / this.region_length;
        if (((this.reconstruction_end-this.reconstruction_start+1) % this.region_length ) !=0) {
        	num_regions=num_regions+1;
        }
        
        this.region_start= new int [num_regions];
        this.region_end= new int [num_regions];
        this.tiling_region_start = new int [num_regions-1 ];
        this.tiling_region_end = new int [num_regions-1 ];
        
        for (int i=0; i< num_regions;i++) {
        	this.region_start[i]= this.reconstruction_start+ i* this.region_length;
        	this.region_end[i]= this.reconstruction_start+ (i+1)* this.region_length-1;
        	if (this.region_end[i]> this.reconstruction_end) {
        		this.region_end[i] = this.reconstruction_end;
        	}
        }
        
        for (int i=0; i< (num_regions-1);i++) {
        	this.tiling_region_start[i]=  (this.region_start[i] + this.region_end[i])/2+1;
        	this.tiling_region_end[i]=  (this.region_start[i+1] + this.region_end[i+1])/2;
        }
        
        BufferedWriter bw = new BufferedWriter(new FileWriter(this.out_dir+"/intermediate/DC_PLAN.txt", false));
        String ss ="Splited REGION:";
        bw.write(ss+"\n");
        ss=Integer.toString(this.region_start[0])+":"+Integer.toString(this.region_end[0]) ;
        for (int i=1; i< num_regions;i++) {
        	ss= ss+"\t"+ Integer.toString(this.region_start[i])+":"+Integer.toString(this.region_end[i]) ;
        }
        bw.write(ss+"\n");
        ss ="TILING_REGION:";
        bw.write(ss+"\n");
        ss=Integer.toString(this.region_start[0])+":"+Integer.toString(this.region_end[0]) ;
        for (int i=1; i< (num_regions-1);i++) {
        	ss= ss+"\t"+ Integer.toString(this.region_start[i])+":"+
        			Integer.toString(this.region_end[i]) ;
        }
        bw.write(ss+"\n");	
        bw.close();
    
	}
}
