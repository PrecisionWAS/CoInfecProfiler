package CoInfecProfiler;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Properties;
import java.util.Scanner;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

public class SARS_CoV_2 {
	public String GATK;
	public String BWA;
	public String samtools;
	public String project_name;
	public String fq_1;
	public String fq_2;
	public String out_dir;
	public String SNP_Matrix;
	
	public  HashMap<String, String> Pos2Ref;
	public  HashMap<String, String> Pos2Alt;
	public  HashMap<String, Double> Pos2AlleleFreq;
	public  HashMap<Integer, Double> Pos2Weights;
	public  HashMap<String, String> Variant2Genotype;
	public  HashMap<String, String> Haplo2Genotype;
	
	public  ArrayList<Integer> PosArr;
	public int  Reconstruction_Start ;
	public int  Reconstruction_End ;
	public double SNV_Cutoff;
	public String [] Potential_Variants;
	public String Regression_Gamma_Min;
	public String Regression_Gamma_Max;
	public String Regression_n_Gamma;
	public String Regression_Lambda;
	public String Reference;
	public double Min_Hap_Freq;
	public double Potential_Variants_Weight;
	public double AF_Zero_Weight;
	public String Rscript_Path;
	public  ArrayList<String> HaploArr;
	public  ArrayList<String> VarArr;
	public String Sum_Weight;
	public int Num_Loci;
	public int Min_Size;
	public int Max_Size = 10;
	public String Regularized_Regression;
	
    public double Regression_MAF_Weight;
    public double Regression_LD_Weight;
    public double	Min_LD_Reads;
    double[][] total_counts;
    double[][] ld_counts;
    public  HashMap<Integer, Integer>  Pos2Index;
    public  HashMap<Integer, Integer>  Index2Pos;
    public double Min_Coverage = 100.0;
    
	public SARS_CoV_2(String parameter_file) throws IOException {
		InputStream is = new FileInputStream(parameter_file);
        Properties prop = new Properties();
        prop.load(is);
		this.BWA = prop.getProperty("BWA");
		this.samtools = prop.getProperty("Samtools");
		this.GATK = prop.getProperty("GATK");
        this.SNP_Matrix = prop.getProperty("SNP_Matrix");
		this.fq_1 = prop.getProperty("Fastq_1");
		this.fq_2 = prop.getProperty("Fastq_2");
        this.out_dir = prop.getProperty("Output_Path")+"/" ;
        this.project_name = prop.getProperty("Proj_Name") ;
        this.Min_Hap_Freq = Double.parseDouble(prop.getProperty("Min_Hap_Freq"));
        this.SNV_Cutoff = Double.parseDouble(prop.getProperty("SNV_Cutoff"));
        this.Reconstruction_Start = Integer.parseInt(prop.getProperty("Reconstruction_Start"));
        this.Reconstruction_End = Integer.parseInt(prop.getProperty("Reconstruction_End"));
        this.Regression_Gamma_Min = prop.getProperty("Regression_Gamma_Min") ;
        this.Regression_Gamma_Max = prop.getProperty("Regression_Gamma_Max") ;
        this.Regression_n_Gamma = prop.getProperty("Regression_n_Gamma") ;
        this.Regression_Lambda = prop.getProperty("Regression_Lambda") ;
        this.Rscript_Path = prop.getProperty("Rscript_Path");
        this.Reference = prop.getProperty("Reference");
        this.Sum_Weight = prop.getProperty("Sum_Weight");
        this.Regularized_Regression = prop.getProperty("Regularized_Regression");
        this.Regression_MAF_Weight = Double.parseDouble(prop.getProperty("MAF_Weight"));
        this.Regression_LD_Weight = Double.parseDouble(prop.getProperty("LD_Weight"));
        this.Min_LD_Reads = Integer.parseInt(prop.getProperty("Min_LD_Reads"));
        
        this.AF_Zero_Weight = Double.parseDouble(prop.getProperty("AF_Zero_Weight"));

        new File(String.valueOf(this.out_dir)).mkdir();
        is.close();
	}

	public void run_bwa() throws IOException, InterruptedException {
		String command = this.BWA+" mem "+this.Reference+" "+this.fq_1+" "+this.fq_2+" > "+this.out_dir+this.project_name+".sam";
		cmd(command);
	}

	public void run_GATK() throws IOException, InterruptedException {
		String command1 = this.samtools+" view -bS "+this.out_dir+this.project_name+".sam"+" > "+this.out_dir+this.project_name+".bam";
		cmd(command1);

		String command2 = this.GATK+" AddOrReplaceReadGroups -I "+this.out_dir+this.project_name+".bam"+" -O "+
				this.out_dir+this.project_name+"_RG.bam"+" -RGID "+this.project_name+" -RGLB "+this.project_name+
				" -RGPL "+this.project_name+" -RGSM "+this.project_name+" -RGPU "+this.project_name;
		cmd(command2);

		String command3 = this.GATK+" SortSam -I "+this.out_dir+this.project_name+"_RG.bam"+" -O "+this.out_dir+this.project_name+"_S.bam -CREATE_INDEX true -VALIDATION_STRINGENCY LENIENT -SO coordinate";
		cmd(command3);

		String command4 = this.GATK+" SplitNCigarReads -R "+this.Reference+" -I "+this.out_dir+this.project_name+"_S.bam -O "+this.out_dir+this.project_name+"_SNCR.bam";
		cmd(command4);

		String command5 = this.GATK+" Mutect2 -R "+this.Reference+" -I "+this.out_dir+this.project_name+"_SNCR.bam -O "+this.out_dir+"known-sites.vcf";
		cmd(command5);

		String command6 = this.GATK+" BaseRecalibrator -R "+this.Reference+" -I "+this.out_dir+this.project_name+"_SNCR.bam --known-sites "+this.out_dir+"known-sites.vcf -O "+this.out_dir+"BQSR.table";
		cmd(command6);

		String command7 = this.GATK+" ApplyBQSR -R "+this.Reference+" -I "+this.out_dir+this.project_name+"_SNCR.bam --bqsr-recal-file "+this.out_dir+"BQSR.table -O "+this.out_dir+this.project_name+"_BQSR.bam";
		cmd(command7);

		String command8 = this.samtools+" index "+this.out_dir+this.project_name+"_BQSR.bam";
		cmd(command8);

		String command9 = this.GATK+" Mutect2 -R "+this.Reference+" -I "+this.out_dir+this.project_name+"_BQSR.bam -mnp-dist 0 -O "+this.out_dir+"mutect2.vcf";
		cmd(command9);
	}
	
	public void readGATKvcf () throws IOException { //Previous GATK VCF format 
		this.Pos2Ref = new HashMap<String, String>  ();
		this.Pos2Alt = new HashMap<String, String>  ();
		this.Pos2AlleleFreq = new HashMap<String, Double> ();
		String line= "";
		BufferedReader br_vcf = new BufferedReader(new FileReader(this.out_dir+"mutect2.vcf"));
		while ((line = br_vcf.readLine()) != null) {
			line= line.replace("\n", "").replace("\r", "");
			if ( (line.length()> 1) &&  (!line.substring(0, 1).equals("#") ) ) {
				String[] tmp = line.split("\t");
				if  (! tmp[4].contains(",") ) {
					this.Pos2Ref.put( tmp[1]+":"+ tmp[3]+";"+tmp[4]
							, tmp[3] );
					this.Pos2Alt.put(tmp[1]+":"+ tmp[3]+";"+tmp[4]
							, tmp[4] );
					String[] tmp2 = tmp[9].split(":");
					String[] tmp3 = tmp2[1].split(",");
					double allele_freq= 0.0; 
					double ref_num = Double.parseDouble(tmp3[0]);
					double alt_num = Double.parseDouble(tmp3[1]);
					
					if ((ref_num+ alt_num) > 1.0 ) {
						allele_freq = (alt_num)/ (ref_num+ alt_num);
						this.Pos2AlleleFreq.put(
								tmp[1]+":"+tmp[3]+";"+ tmp[4] , allele_freq);
					}
					
				} else {
					String[] tmp4 = tmp[4].split(","); //multiple alt alleles, such as A,T *,G
					String[] tmp5 = tmp[9].split(":"); 
					String[] tmp6 = tmp5[1].split(","); // G	A,T counts: 100,200,3000
					
					double total_count=0.0;
					
					for (int j=0;j < tmp6.length; j++) {
						total_count = total_count + 
								Double.parseDouble(tmp6[j]) ;
					}
					
					if ((total_count) > 1.0 ) {
						for (int j=1;j < tmp6.length   ; j++) {
							this.Pos2Ref.put( tmp[1]+":"+ tmp[3]+";"+tmp[4] , tmp[3] );
							this.Pos2Alt.put( tmp[1]+":"+ tmp[3]+";"+tmp[4], tmp[4]  );
							
							this.Pos2AlleleFreq.put( tmp[1]+ ":" +tmp[3]+";"+ tmp4[j-1]   , 
									Double.parseDouble(tmp6[j]  )/ total_count);					
						}
					}
					
				}
			}
		}
		br_vcf.close();
	}

	public void readvcf() throws IOException { //MUTECT2 VCF format
		this.Pos2Ref = new HashMap<String, String>  ();
		this.Pos2Alt = new HashMap<String, String>  ();
		this.Pos2AlleleFreq = new HashMap<String, Double> ();
		String line= "";
		int COUNT=0;
		BufferedReader br_vcf = new BufferedReader(new FileReader(this.out_dir+"mutect2.vcf"));
		while ((line = br_vcf.readLine()) != null) {
			line= line.replace("\n", "").replace("\r", "");
			String[] tmp = line.split("\t");
			boolean MinCovPASS=false;
			if ( (line.length()> 1) &&  (!line.substring(0, 1).equals("#") ) ) {
				String[] tmp55 = tmp[9].split(":"); 
				String[] tmp66 = tmp55[1].split(","); // G	A,T counts: 100,200,3000
				
				double total_count=0.0;
				
				for(int j=0;j < tmp66.length; j++) {
					total_count = total_count + 
							Double.parseDouble(tmp66[j]) ;
				}
				if (total_count> this.Min_Coverage) {
					MinCovPASS=true;
				}
			}
		
				 
			
			if ( ( (line.length()> 1) &&  (!line.substring(0, 1).equals("#") ) )  && 
			MinCovPASS ){
				if  (! tmp[4].contains(",") ) {				
						this.Pos2Ref.put( tmp[1]+":"+ tmp[3]+";"+tmp[4] , tmp[3] );
						this.Pos2Alt.put( tmp[1]+":"+ tmp[3]+";"+tmp[4] , tmp[4] );
						String[] tmp2 = tmp[9].split(":");
						String[] tmp3 = tmp2[2].split(",");
						double allele_freq=Double.parseDouble(tmp3[0]);  
						this.Pos2AlleleFreq.put(tmp[1]+":"+ tmp[3]+";"+ tmp[4] , allele_freq);
					
					
				} else {
					if ( tmp[3].length()==1 ) {
						String[] tmp4 = tmp[4].split(","); //multiple alt alleles, such as A,T *,G
						String[] tmp5 = tmp[9].split(":"); 
						String[] tmp6 = tmp5[2].split(","); // G	A,T counts: 100,200,3000
						
						if (tmp6.length==tmp4.length){
//							System.out.println(line);
							for (int j=0;j < tmp4.length  ; j++) {
//								System.out.println(this.removeredundance(tmp[3], tmp4[j]));
									
								
								this.Pos2Ref.put( tmp[1]+":"+ tmp[3]+";"+tmp[4] , tmp[3] );
								this.Pos2Alt.put( tmp[1]+":"+ tmp[3]+";"+tmp[4], tmp4[j]  );
								this.Pos2AlleleFreq.put( tmp[1]+ ":" +tmp[3]+";"+ tmp4[j]   , 
											Double.parseDouble(tmp6[j]  ));		
//								System.out.println(tmp6[j] );		
							}
						}else {
							System.out.println("Error:\t"+line);
						}
					}else{
						String[] tmp4 = tmp[4].split(","); //multiple alt alleles, such as A,T *,G
						String[] tmp5 = tmp[9].split(":"); 
						String[] tmp6 = tmp5[2].split(","); // G	A,T counts: 100,200,3000
						if (tmp6.length==tmp4.length){
							int initPOS= Integer.parseInt(tmp[1]);
//							System.out.println(line);
							for (int j=0;j < tmp4.length  ; j++) {
								if (tmp4[j].length() == tmp[3].length()) {
//									System.out.println(line);
									for (int p=0;p< tmp4[j].length();p++) {
										if (!tmp4[j].substring(p,p+1).
												equals(tmp[3].substring(p,p+1) ) ) {
											String key= Integer.toString(initPOS+p)+":"+tmp[3].substring(p,p+1)
													+";"+tmp4[j].substring(p,p+1);
											this.Pos2Ref.put( key, tmp[3].substring(p,p+1) );
											this.Pos2Alt.put( key, tmp4[j].substring(p,p+1) );
											this.Pos2AlleleFreq.put( Integer.toString(initPOS+p )+ ":" 
													+tmp[3].substring(p,p+1)+";"+ tmp4[j].substring(p,p+1) , 
													Double.parseDouble(tmp6[j]  ));	
//											System.out.println(Integer.toString(initPOS+p )+"\t"+
//													 tmp[3].substring(p,p+1)+"\t"+
//													 tmp4[j].substring(p,p+1)+"\t"+
//													 tmp6[j]);
//											if (p>0) {
//												System.out.println(initPOS+p);	
//												System.out.println(tmp[3].substring(p,p+1));
//											}
										}
									}
								}else {
//									System.out.println(line);
									if (  (tmp[3].length()==1) || (tmp4[j].length()==1)) {
//										System.out.println(line);
										this.Pos2Ref.put( tmp[1]+ ":" +tmp[3]+";"+tmp4[j], tmp[3] );
										this.Pos2Alt.put(tmp[1]+ ":" +tmp[3]+";"+tmp4[j] , tmp4[j]  );
										this.Pos2AlleleFreq.put( tmp[1]+ ":" +tmp[3]+";"+tmp4[j]   , 
													Double.parseDouble(tmp6[j]  ));	
//										System.out.println(tmp[1]+"\t"+tmp[3]+"\t"+tmp4[j]+"\t"
//												+tmp[1]+ ":" +tmp4[j]+"\t"+tmp6[j] );
									} else {
										int remove =removeredundance (tmp[3], tmp4[j]);
//										if (remove>1) {
//											remove=1;
//										}
										
										String key= Integer.toString(initPOS+ remove)+ ":"+
												tmp[3].substring(remove) +";"+
												tmp4[j].substring(remove);
										this.Pos2Ref.put( key, tmp[3].substring(remove) );
										this.Pos2Alt.put( key, tmp4[j].substring(remove) );
										this.Pos2AlleleFreq.put( key , 
												Double.parseDouble(tmp6[j]  ));	
										
//										System.out.println(line);
//										System.out.println(key+"\t"+ tmp6[j] );
										COUNT++;
										
									}
								}
							}
						}
					}
				}
			}
		}
		br_vcf.close();
//		System.out.println(COUNT);
	}

	public void variantmaker () throws IOException {
		int count=0;
		BufferedWriter bw = new BufferedWriter(new FileWriter(
				this.out_dir+"/snv.txt", false));
		for (Object key:this.Pos2Ref.keySet()) {
			String var = key.toString();
			if ((var.split(":")[1].split(";")[0].length()==1)  
					&& (var.split(":")[1].split(";")[1].length()==1) ){
				count++;
				bw.write(var.split(":")[0] +"\t"+ var.split(":")[1].split(";")[0]+
						"\t"+var.split(":")[1].split(";")[1] +"\n");
			}
		}
//		for (int i= this.Reconstruction_Start; i<= this.Reconstruction_End;i++) {
//			if  (this.Pos2Ref.containsKey(i) ){
//				if ((this.Pos2Ref.get(i).length()==1)  && (this.Pos2Alt.get(i).length()==1) ){
//					count++;
//					bw.write(Integer.toString(i) +"\t"+ this.Pos2Ref.get(i)+"\t"+this.Pos2Alt.get(i) +"\n");
//				}
//			}
//		}

		bw.close();
		this.total_counts = new double [count][count];
		this.ld_counts = new double [count][count];
		this.Pos2Index = new HashMap<Integer, Integer>  ();
		this.Index2Pos = new HashMap<Integer, Integer>  ();
		int index =0;
		for (int i= this.Reconstruction_Start; i<= this.Reconstruction_End;i++) {
			for (Object key:this.Pos2Ref.keySet()) {
				String var = key.toString();
//				System.out.println(var.split(":")[0]+"\t"+Integer.toString(pos1));
				if ( Integer.parseInt((var.split(":")[0]) )==i ) {
					if( var.split(":")[1].length()==3) {
						this.Pos2Index.put(i, index);
						this.Index2Pos.put(index, i);
						index++;
					}
				}
			}
//			if  (this.Pos2Ref.containsKey(i) ){
//				if ((this.Pos2Ref.get(i).length()==1)  && (this.Pos2Alt.get(i).length()==1) ){
//					this.Pos2Index.put(i, index);
//					this.Index2Pos.put(index, i);
//					index++;
//				}
//			}
		}
		for (int i=0;i < this.Pos2Index.size();i++) {
			for (int j =0;j < this.Pos2Index.size();j++) {
				this.total_counts[i][j]=0;
				this.ld_counts[i][j]=0;
			}
		}
		
		
	}
	
	public void vefgenerater () throws IOException {
		HashSet<Integer> trueVarPos = new HashSet<Integer>(); 
		HashMap<Integer,HashMap<String,VarObj>> variantEncoder = 
				new HashMap<Integer,HashMap<String,VarObj>>(); 
		String line="";
		int variantCode=0;
		BufferedReader br_snv = new BufferedReader
				(new FileReader(this.out_dir+"/snv.txt"));
		while ((line = br_snv.readLine()) != null) {
				line= line.replace("\n", "").replace("\r", "").replace(" ", "");
				String[] tmp = line.split("\t");
				if  (!tmp[2].equals("*") ) {
					int  pos = Integer.parseInt(tmp[0]);
					String ref = tmp[1];
//					String alt =  tmp[2];
					trueVarPos.add(Integer.parseInt(tmp[0])); 	
					variantEncoder.put(pos, new HashMap<String,VarObj>());
					VarObj altAlleleAtPos = new VarObj(pos,variantCode);
					variantEncoder.get(pos).put(ref,altAlleleAtPos);
				}
		}
		br_snv.close();
		VEFMaker(this.out_dir+this.project_name+".sam", this.out_dir+"/raw.vef", variantEncoder);
		PairedReadLinker.link_paired_vef(this.out_dir+"/raw.vef", this.out_dir+"/filter.vef");
	}
	
	public int getposition (String s1) throws IOException {
		String[] tmp = s1.split("=");
		return Integer.parseInt(tmp[0]);
	}
	
	public int removeredundance(String s1, String s2) throws IOException {
		int minlen=s1.length();
		if (s2.length()< minlen) {
			minlen= s2.length();
		}
		int pos=0;
		for (int i=0;i<(minlen -1);i++) {
			if ( s1.substring(i, i+1).equals(s2.substring(i, i+1) )) {
				pos++;
			} else {
				break;
			}
		}
		return pos;
		
	}
	
	public int getallele (String s1) throws IOException {
		String[] tmp = s1.split("=");
		return Integer.parseInt(tmp[1]);
	}
	
	public void readref () throws IOException {
		BufferedReader br_vef = new BufferedReader
				(new FileReader(this.out_dir+"/filter.vef"));
		String line="";
		while ((line = br_vef.readLine()) != null) {
			
			line= line.replace("\n", "").replace("\r", "").replace(" ", "");
			String[] tmp = line.split("\t");
			String[] tmp2 = tmp[1].split(";");
			if (tmp2.length>1) {
				for (int i=0;i< tmp2.length;i++) {
					for (int j= (i+1);j< tmp2.length;j++) {
//						System.out.println( line);
						int pos1= getposition (tmp2[i]);
						int pos2= getposition (tmp2[j]);
						int allele1 =getallele (tmp2[i]);
						int allele2 =getallele (tmp2[j]);
						if (this.Pos2Index.containsKey(pos1) && this.Pos2Index.containsKey(pos2)) {
							int index1= this.Pos2Index.get(pos1);
							int index2= this.Pos2Index.get(pos2);
							this.total_counts[index1][index2] = 
									this.total_counts[index1][index2]+ 1.0;
							if  ( (allele1==1 ) && (allele2==1 )) {
								this.ld_counts[index1][index2] = 
										this.ld_counts[index1][index2]+ 1.0;
							}
						}
						
					}
				}
			}
		}
		
		br_vef.close();
		
		
	}

	public static void VEFMaker(String input_SAM, String output_vef_raw,
			HashMap<Integer,HashMap<String,VarObj>> variantEncoder ) 
			throws FileNotFoundException {
		
		String seq_tech ="linked_reads";
		Set<Integer> keyMATCH, keyINS, keyDEL, variantKS = variantEncoder.keySet(), tempSet, indelKS;
		SortedSet<Integer> variantSS = new TreeSet<Integer>(); 
		variantSS.addAll(variantKS); 
		File BAMPath = new File(input_SAM); 
		Scanner BAMScanner = new Scanner(BAMPath);  
		

		BAMScanner.useDelimiter("\t");
		String currLine = BAMScanner.nextLine(); 
		// System.out.println(currLine);
		while (currLine.matches("@(.*)")) {
			// String prevLine = currLine;
			currLine = BAMScanner.nextLine(); 
			// if (prevLine.contains("PN:")) break;
			// System.out.println(currLine);
		}
		String QNAME, CIGAR, SEQ, tempPos = "",currAltAllele;
		Integer FLAG, startPOS = 1, currentPos, posToAdd, matchToIndel, skipSBases = 0;
		VarObj currVarObj; 
		HashMap<Integer, HashMap<Integer,String>> hmMATCH = new HashMap<Integer, HashMap<Integer,String>>();
		HashMap<Integer, HashMap<Integer,String>> hmINS = new HashMap<Integer, HashMap<Integer,String>>();
		HashMap<Integer, HashMap<Integer,String>> hmDEL = new HashMap<Integer, HashMap<Integer,String>>();
		HashMap<Integer,String> defaultHM, tempFinderHM;
		HashMap<String,VarObj> tempReporterHM; 
		/* HashMap<Integer,Integer> hmVarCount = new HashMap<Integer,Integer>(); 
		for (Integer i : variantSS) {
			hmVarCount.put(i,0);
		}
		int readCount = 0; */
		//PrintWriter VEFFile = new PrintWriter(BAMPrefix + ".raw.vef");// replaced with the line below
		PrintWriter VEFFile = new PrintWriter(output_vef_raw);
		while (BAMScanner.hasNextLine()) {
			// readCount++; 
			QNAME = BAMScanner.next();
			String QNAME_tmp = QNAME;
			if (seq_tech.equals("10x_linked-reads")) {
//				System.out.println(seq_tech);
				QNAME="";
			}
			// System.out.println(QNAME);
			FLAG = BAMScanner.nextInt();
			if (FLAG > 163) {
				BAMScanner.nextLine();
				continue; 	// If the read is chimeric, or mapped to more than one location, or of poor technical quality, skip it. 
			}
			BAMScanner.next(); 	// Skip the FLAG, RNAME.
			startPOS = BAMScanner.nextInt();
			currentPos = startPOS; 
			BAMScanner.next(); 		// Skip the MAPQ.			
			CIGAR = BAMScanner.next();
			if (CIGAR.charAt(0) == '*') {
				BAMScanner.nextLine(); 	// * refers to the fact that no CIGAR information is available. 
				continue;	 			// The above line skips the rest of the information in this read entirely.
			}
			for (int c = 0; c < CIGAR.length(); c++) {
				Character tempChar = CIGAR.charAt(c); 
				if (tempChar.compareTo('M') == 0) {
					posToAdd = Integer.parseInt(tempPos);
					for (int addPos = currentPos; addPos < currentPos + posToAdd; addPos++) {
						defaultHM = new HashMap<Integer,String>(); 
						defaultHM.put(1, "N"); 
						hmMATCH.put(addPos,defaultHM);
					}
					currentPos += posToAdd;
					tempPos = ""; 
				} else if (tempChar.compareTo('I') == 0) { 
					// Command to check for insertions: samtools view HIV_p1_sample_1.procd.bam | awk '{if ($6 ~ /I/) print $6;}' -
					posToAdd = Integer.parseInt(tempPos);
					matchToIndel = currentPos - 1; 	// This is the first position of the insertion, which it is recognized by. 
					hmMATCH.remove(matchToIndel); 	// Remove the first position of the insertion from the match HM.
					defaultHM = new HashMap<Integer,String>(); 
					defaultHM.put(posToAdd + 1, "N");	// This accounts for the entire length of the insertion.
					hmINS.put(matchToIndel, defaultHM); 
					currentPos = matchToIndel + 1;		// Since the insertion starts at the previous currentPos - 1, 
				//the next match along starts at matchToIndel + 1 = previous currentPos.
					tempPos = ""; 
				} else if (tempChar.compareTo('D') == 0) {
					posToAdd = Integer.parseInt(tempPos);
					matchToIndel = currentPos - 1; 	// This is the first position of the deletion, which it is recognized by. 
					hmMATCH.remove(matchToIndel); 	// Remove the first position of the deletion from the match HM.
					defaultHM = new HashMap<Integer,String>(); 
					defaultHM.put(posToAdd + 1, "N");	// This accounts for the entire length of the deletion.
					hmDEL.put(matchToIndel, defaultHM); 
					currentPos += posToAdd; 	// Since the deletion starts at the previous currentPos - 1, the next match along starts at currentPos + matchToIndel - 1.
					tempPos = ""; 
				} else if (tempChar.compareTo('S') == 0) {
					skipSBases = Integer.parseInt(tempPos);	// These bases do appear in SEQ BUT they are not aligned to the reference genome. Therefore these bases do not contribute to the currentPOS. 
					tempPos = ""; 				
				} else if (tempChar.compareTo('H') == 0) {
					tempPos = ""; // These bases DO NOT appear in SEQ and they are not aligned to the reference genome. Therefore these bases do not contribute to the currentPOS. 				
				} else { 
					tempPos += tempChar;
				}
			}
			
			keyMATCH = hmMATCH.keySet(); 
			keyINS = hmINS.keySet();
			keyDEL = hmDEL.keySet(); 
			for (int s = 0; s < 3; s++) {
				BAMScanner.next(); 	// Skip the RNEXT, PNEXT, TLEN.
			}
			SEQ = BAMScanner.next();
			int endPOS = startPOS + SEQ.length() - 1;
			currentPos = startPOS; 
			// 4. Test the following to see if the SEQ-parsing code works properly, and the proper bases are assigned to the correct Hash Maps. 
			for (int s = skipSBases; s < SEQ.length(); s++) {	// Skip all of the bases that were soft-clipped but still reported at the beginning of SEQ. In the event there is no 'S' in the CIGAR string, this is 0 and we start reporting from the end of SEQ.  
				// System.out.println(s + "\t" + currentPos);
				if (keyMATCH.contains(currentPos)) {
					hmMATCH.get(currentPos).put(1,SEQ.substring(s,s+1));
					currentPos++; 
				} else if (keyINS.contains(currentPos)) {
					tempSet = hmINS.get(currentPos).keySet();
					for (Integer length : tempSet) {
						hmINS.get(currentPos).put(length,SEQ.substring(s,s+length)); 
					}
					currentPos++; 
				} else if (keyDEL.contains(currentPos)) {
					tempSet = hmDEL.get(currentPos).keySet();
					for (Integer length : tempSet) {
						hmDEL.get(currentPos).put(length,SEQ.substring(s,s+1));
						currentPos += length; 
					}
				}
			}

			skipSBases = 0; // Reset the 'start' position of SEQ-processing to 0 in case there aren't any S-es in the next CIGAR string. 
			
			StringBuilder readInfo = new StringBuilder(); 
			// System.out.println(startPOS + "\t" + endPOS);
			for (Integer VCFVarPos : variantSS.subSet(startPOS,endPOS+1)) {	// subSet method is (inclusive,exclusive)
				// System.out.println(VCFVarPos); 
				tempReporterHM = variantEncoder.get(VCFVarPos);
				if (keyMATCH.contains(VCFVarPos)) {
					currAltAllele = hmMATCH.get(VCFVarPos).get(1);
					// System.out.println(currAltAllele);
					if (tempReporterHM.containsKey(currAltAllele)) {
						currVarObj = tempReporterHM.get(currAltAllele);
						// VEFFile.print(currVarObj.intM + "=" + currVarObj.varCode + ";");
						readInfo.append(currVarObj.intM + "=" + currVarObj.varCode + ";");
						// hmVarCount.put(VCFVarPos, hmVarCount.get(VCFVarPos) + 1);
					} else {
						readInfo.append(VCFVarPos + "=1;"); // Fixed as of 10282018. Noted by CC that in some reads, variant positions in those reads were not being annotated.
						//System.out.println(VCFVarPos + " done");	// Realized that the 'else' clause was in the wrong place i.e.: alleles  =/= alternate (i.e: the reference) 
					}												// were not noted when they existed. 
					// System.out.println("match");
				} else if (keyINS.contains(VCFVarPos)) {
					tempFinderHM = hmINS.get(VCFVarPos);
					indelKS = tempFinderHM.keySet(); 
					for (Integer indel : indelKS) {
						currAltAllele = tempFinderHM.get(indel);
						if (tempReporterHM.containsKey(currAltAllele)) {
							currVarObj = tempReporterHM.get(currAltAllele);
							// VEFFile.print(currVarObj.intM + "=" + currVarObj.varCode + ";");
							readInfo.append(currVarObj.intM + "=" + currVarObj.varCode + ";");
							// hmVarCount.put(VCFVarPos, hmVarCount.get(VCFVarPos) + 1);
						} else {
							readInfo.append(VCFVarPos + "=1;");
							// System.out.println(VCFVarPos + " done");
						}
					}
					// System.out.println("ins");
				} else if (keyDEL.contains(VCFVarPos)) {
					tempFinderHM = hmDEL.get(VCFVarPos);
					indelKS = tempFinderHM.keySet(); 
					for (Integer indel : indelKS) {
						currAltAllele = tempFinderHM.get(indel);
						if (tempReporterHM.containsKey(currAltAllele)) {
							currVarObj = tempReporterHM.get(currAltAllele);
							// VEFFile.print(currVarObj.intM + "=" + currVarObj.varCode + ";");
							readInfo.append(currVarObj.intM + "=" + currVarObj.varCode + ";");
							// hmVarCount.put(VCFVarPos, hmVarCount.get(VCFVarPos) + 1);
						} else {
							readInfo.append(VCFVarPos + "=1;");
							// System.out.println(VCFVarPos + " done");
						}
						// System.out.println("del");
					}
				}
			}
			// System.out.println(readInfo.toString());
			if (readInfo.toString().isEmpty()) {
				BAMScanner.nextLine(); 	// 
				hmMATCH.clear();
				hmINS.clear();
				hmDEL.clear(); 
				continue;
			}
			// System.out.println(readInfo);
			int flag=0;
			if (seq_tech.equals("10x_linked-reads")) {
				String DM= "";
//				OM:i:60 XM:Z:0  TQ:Z:>@>>=?@    TR:Z:TAGGGTT    AS:i:0  XS:i:-70        XT:i:0  BX:Z:TGTCACCAGGGTACGT-1 
				String tmp ="";
				String BX="";
				for (int s = 0; s < 14; s++) {
					tmp = BAMScanner.next(); 	
					if ((tmp.length()>4) && (tmp.substring(0, 3).equals("BX:")) ) {
						BX= tmp;
						flag++;
					}
					if ((tmp.length()>4) && (tmp.substring(0, 3).equals("DM:")) ) {
						DM= tmp;
						flag++;
					}
					if (flag==2) {
						break;
					}
				}
				QNAME=DM.replace("\n", "")+":"+ BX.replace("\n", "");
			}
			if (flag==2 ) {
				String[] x_arr =  QNAME_tmp.split("_");
				String hap = x_arr[0]+"_"+ x_arr[1]+ "_"+ x_arr[2];
				VEFFile.println(QNAME+"_"+hap + ":\t" + readInfo +"\t"+QNAME_tmp+ "\t"+ startPOS+"\t" + endPOS);
			} 
			if (!seq_tech.equals("10x_linked-reads")) {
				VEFFile.println(QNAME + ":\t" + readInfo +"\t"+QNAME_tmp+ "\t"+ startPOS+"\t" + endPOS);
			}
			
			BAMScanner.nextLine(); 
			hmMATCH.clear();
			hmINS.clear();
			hmDEL.clear(); 
		}
		VEFFile.close();
		BAMScanner.close();
//		System.out.println(output_vef_raw +" has been generated.");
	}
	
	public void readmatrix () throws IOException {
		this.Variant2Genotype = new HashMap<String, String> ();
		HashMap<Integer, String> index2variant = new HashMap<Integer, String>();
		this.Haplo2Genotype = new HashMap<String, String> ();
		this.HaploArr = new ArrayList<String>();
		this.VarArr = new ArrayList<String>();
		String line= "";
		BufferedReader br_matrix = new BufferedReader(new FileReader(this.SNP_Matrix));
		while ((line = br_matrix.readLine()) != null) {
			line= line.replace("\n", "").replace("\r", "");
			if ( (line.length()> 1) &&  (!line.substring(0, 1).equals("#") ) ) {
				String[] tmp = line.split("\t");
				if (!line.substring(0, 3).equals("SNP") ) {
					this.HaploArr.add(tmp[0]) ;
					this.Haplo2Genotype.put(tmp[0], "");
					for (int i=1; i< tmp.length; i++) {
						String var = index2variant.get(i);
						this.Variant2Genotype.put(var, this.Variant2Genotype.get(var)+ tmp[i] );
						this.Haplo2Genotype.put(tmp[0] ,this.Haplo2Genotype.get(tmp[0])+ tmp[i] );
					}
//					for (int i=0;i < this.Potential_Variants.length;i++) {
//						if (tmp[0].equals(this.Potential_Variants[i]){
//					}
				}
				if (line.substring(0, 3).equals("SNP") ) {
					for (int i=1; i< tmp.length; i++) {
						this.Variant2Genotype.put(tmp[i], "");
						index2variant.put(i, tmp[i]);
						this.VarArr.add(tmp[i]);
					}	
				}
			}
		}
//		System.out.println(this.Variant2Genotype.get("21:C|T") );
		br_matrix.close();
	}
	
	public void l0l1init () throws IOException {
		int loci_num= 0;
		new File(String.valueOf(this.out_dir) ).mkdir();
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(
				this.out_dir+"/matrix_init.txt", false));
		String ss="0|";
		for (int i=0 ;i <= this.HaploArr.size();i++ ) {
			ss = ss+" "+ this.Sum_Weight;
		}
		bw.write(ss+"\n");
		for (Object key:this.Pos2Ref.keySet()) {
			String var = key.toString();
			
			if ( this.Variant2Genotype.containsKey(var) ) {
//				System.out.println(var);
				ss=var +"| "+ Double.toString(this.Pos2AlleleFreq.get(var)*this.Regression_MAF_Weight);
				for (int j=0; j< this.Variant2Genotype.get(var).length(); j++) {
					ss=ss+" "+Double.toString (Double.parseDouble(this.Variant2Genotype.get(var).
							substring(j, j+1))*this.Regression_MAF_Weight);
				}
				loci_num++;
				bw.write(ss+"\n");
			}
			
		}
		
		for (int i= this.Reconstruction_Start; i<= (this.Reconstruction_End*-1 );i++) {// delete
			
			if (this.Pos2Ref.containsKey(i)) {
				String var = Integer.toString (i)+":"+ 
						this.Pos2Ref.get(i)+"|" + this.Pos2Alt.get(i);
				
				if  (! this.Pos2Alt.get(i).contains(",") ) {
					if ( this.Variant2Genotype.containsKey(var) ) {
						ss=var +"| "+ Double.toString(this.Pos2AlleleFreq.get(Integer.toString (i)+":"+ 
								this.Pos2Ref.get(i)+":"+this.Pos2Alt.get(i))*this.Regression_MAF_Weight);
						for (int j=0; j< this.Variant2Genotype.get(var).length(); j++) {
							ss=ss+" "+Double.toString (Double.parseDouble(this.Variant2Genotype.get(var).
									substring(j, j+1))*this.Regression_MAF_Weight);
						}
						loci_num++;
						bw.write(ss+"\n");
					}
				} else {
					String[] tmp = this.Pos2Alt.get(i).split(",");
					for (int j=0; j < tmp.length; j++) {
						if ( ! tmp[j].equals("*") ) {
							var = Integer.toString (i)+":"+ 
									this.Pos2Ref.get(i)+"|" + tmp[j];
							if ( this.Variant2Genotype.containsKey(var) ) {
								ss=var +"| "+ Double.toString(this.Pos2AlleleFreq.get(Integer.toString (i)+":"+ 
										this.Pos2Ref.get(i)+ ":"+tmp[j])*this.Regression_MAF_Weight);
								for (int k=0; k< this.Variant2Genotype.get(var).length(); k++) {
									ss=ss+" "+Double.toString (Double.parseDouble(this.Variant2Genotype.get(var).
											substring(k, k+1))*this.Regression_MAF_Weight);
								}
								loci_num++;
								bw.write(ss+"\n");
							}
						} else {
							String calgary= "calgary";
						}
					}
				}
			}
		}

		
		
		
		for (int p=0; p< this.VarArr.size();p++) {
			boolean flag=true;
			String var = this.VarArr.get(p);
			String tmp []= this.VarArr.get(p).split(":");
			int pos = Integer.parseInt(tmp[0]);
			for (Object key:this.Pos2Ref.keySet()) {
				String vcfvar = key.toString();
				if (var.equals(vcfvar) ) {
					flag=false;
				}
			}
			if (flag) {
				ss=var +"| "+ Double.toString(0.0*this.AF_Zero_Weight);
				for (int j =0; j< this.Variant2Genotype.get(var).length(); j++) {
					ss=ss+" "+Double.toString (Double.parseDouble(this.Variant2Genotype.get(var).
							substring(j, j+1))*this.AF_Zero_Weight);
				}
				loci_num++;
				bw.write(ss+"\n");
			}
		}
		
			
		
		for (int i= this.Reconstruction_Start; i<= (-1*this.Reconstruction_End );i++) {
			for (int p=0; p< this.VarArr.size();p++) {
				
				String var = this.VarArr.get(p);
				String tmp []= this.VarArr.get(p).split(":");
				int pos = Integer.parseInt(tmp[0]);
				boolean flag=true;
				if  (!  this.Pos2Ref.containsKey(i)) {
					flag=false;
				} else {

					if (! (tmp[1].split("\\|")[0].equals(this.Pos2Ref.get(i)))) {
						flag=false;
					}
					if (!(tmp[1].split("\\|")[1].equals(this.Pos2Alt.get(i)))) {
						flag=false;
					}	
				}
				
				if  ( (i == pos)  &&  (this.Variant2Genotype.containsKey(var))  &&
						( !flag )  ) {
					
					ss=var +"| "+ Double.toString(0.0*this.AF_Zero_Weight);
					for (int j =0; j< this.Variant2Genotype.get(var).length(); j++) {
						ss=ss+" "+Double.toString (Double.parseDouble(this.Variant2Genotype.get(var).
								substring(j, j+1))*this.AF_Zero_Weight);
					}
					loci_num++;
					bw.write(ss+"\n");
				}	
			}
		}
		
		
		
//		for (int i=0;i<this.Pos2Index.size();i++) {
		for (Integer key1:this.Pos2Index.keySet() ) {
//			for (int j= (i+1);j <this.Pos2Index.size();j++) {
			for (Integer key2:this.Pos2Index.keySet() ) {
				int posi= key1.intValue();
				int posj= key2.intValue();
				if (posi==posj) {
					break;
				}
//				System.out.println( Pos2Index.get(posi) );
				int i= this.Pos2Index.get(posi);
				int j= this.Pos2Index.get(posj);
				if ((this.total_counts[i][j]> this.Min_LD_Reads ) && 
						(this.total_counts[i][j]> 1.0 ) ) {
					if ( this.ld_counts[i][j]> 5.0) {
//						System.out.println( Pos2Index.get(i) );
						if ( (this.ld_counts[i][j]/
						this.total_counts[i][j] ) > this.Min_Hap_Freq ) {
							ss=Integer.toString(this.Index2Pos.get(i))+":"+ 
									Integer.toString(this.Index2Pos.get(j))+"| "+
									Double.toString (this.Regression_LD_Weight* this.ld_counts[i][j]/
											this.total_counts[i][j] );
							int pos1= this.Index2Pos.get(i);
							int pos2= this.Index2Pos.get(j);
							String var1="";
							String var2="";
							int  flag=0;
							
							for (Object key:this.Pos2Ref.keySet()) {
								String var = key.toString();
//								System.out.println(var.split(":")[0]+"\t"+Integer.toString(pos1));
								if ( Integer.parseInt((var.split(":")[0]) )==pos1) {
									if( var.split(":")[1].length()==3) {
										flag++;
										var1= var;
										break;
									}
								}
							}
							
							for (Object key:this.Pos2Ref.keySet()) {
								String var = key.toString();
								if ( Integer.parseInt((var.split(":")[0]) )==pos2) {
									if( var.split(":")[1].length()==3) {
										flag++;
										var2= var;
										break;
									}
								}
							}
							
							
							
//							String var1 = Integer.toString (pos1)+":"+ 
//									this.Pos2Ref.get(pos1)+";" + this.Pos2Alt.get(pos1);
//							String var2 = Integer.toString (pos2)+":"+ 
//									this.Pos2Ref.get(pos2)+";" + this.Pos2Alt.get(pos2);
							if (( this.Variant2Genotype.containsKey(var1) && 
									this.Variant2Genotype.containsKey(var2)  ) && (flag==2) ) {
								for (int k=0; k< this.Variant2Genotype.get(var1).length(); k++) {
									ss=ss+" "+Double.toString (Double.parseDouble(this.Variant2Genotype.get(var1).
											substring(k, k+1))*this.Regression_LD_Weight *
											Double.parseDouble(this.Variant2Genotype.get(var2).
													substring(k, k+1)) );
								}
								bw.write(ss+"\n");
							}
						}
						
					}
				}
			}
		}
		
		bw.close();
		this.Num_Loci = loci_num;
	}
	
	public void rscript() throws IOException, InterruptedException {
		String rfile= this.out_dir+ this.project_name+ ".R";
		BufferedWriter bw = new BufferedWriter(new FileWriter( rfile, false));
		bw.write("library(\'L0Learn\')\n");
		bw.write("maxSize= "+Integer.toString(this.Max_Size)+"\n");
		bw.write("minSize= "+Integer.toString(this.Min_Size)+"\n");
		bw.write("mum_sites="+Integer.toString(this.Num_Loci)+"\n");
		bw.write("maf_weights= "+Double.toString( 1.0)+"\n");
		String tmp ="txt<-as.matrix(read.table(\'";
		String infile= this.out_dir+"/matrix_init.txt";
		String outfile= this.out_dir+"/matrix_out.txt";
		
		tmp=tmp+infile+"\',sep=\' \'))";
		bw.write(tmp+"\n");
		bw.write("y<-as.numeric(txt[,2])\n");//apply(X,2,as.numeric)
		bw.write("X<-apply(txt[,3:ncol(txt)], 2,as.numeric)  \n");
		bw.write("if ((ncol(X))>minSize){\n");
		tmp= "cvfit = L0Learn.cvfit(X, y, nFolds=5, seed=2, penalty=\""+this.Regularized_Regression+"\", "
				+ "nGamma="+this.Regression_n_Gamma+", gammaMin="+ 
				this.Regression_Gamma_Min+","
				+ " gammaMax="+ this.Regression_Gamma_Max+ 
				", maxSuppSize=maxSize"+", intercept= FALSE, algorithm=\"CDPSI\")";
		bw.write(tmp+"\n");
			
		bw.write("optimalLambda  = "+ this.Regression_Lambda+ "\n");
		bw.write("optimalGammaIndex= which(unlist(lapply(cvfit$cvMeans, min)) == "
		+ "min(unlist(lapply(cvfit$cvMeans, min))))\n");
		bw.write("tmp =coef(cvfit, lambda=optimalLambda, gamma=cvfit$fit$gamma[optimalGammaIndex])\n");	
		
		tmp= "write( paste( \"#Hap_ID\", \"Freq\" ,\"Haplotype\",sep = \"\\t\" ) "
				+ ",file=\""+outfile+"\",append=FALSE)";
		bw.write(tmp+"\n");
		bw.write("for(i in 1:length(as.vector(tmp@x))){\n");
		
		bw.write("hap = gsub(\", \",\"\",toString(round(X[, "
				+ "as.vector(tmp@i)[i]+1][2:(1+ mum_sites)]/maf_weights)))\n");
		tmp= "write(paste( paste( \"h\", toString(as.vector(tmp@i)[i]),sep = \"\" ), "
				+ "as.vector(tmp@x)[i],hap, sep = \"\\t\"),file=\""+outfile
				+"\",append=TRUE)";
		bw.write(tmp+"\n");
		bw.write("}}\n");
		bw.write("print (\"FINISH\")\n");
		
		bw.close();
	}
	
	public void run_R() throws IOException, InterruptedException {
		ProcessBuilder CMDLine = new ProcessBuilder(this.Rscript_Path, 
				this.out_dir+ "/"+this.project_name+ ".R" );
        Process CMDProcess = CMDLine.start();
            
		BufferedReader br = new BufferedReader
				(new InputStreamReader(CMDProcess.getInputStream()));
		String line;
		while ((line = br.readLine()) != null) {
			line= line.replace("\n", "").replace("\r", "");
		}
		CMDProcess.waitFor();
		br.close();
	}

	public void userformat() throws IOException {
		String infil = this.out_dir+"/matrix_out.txt";
		String hapfil = this.out_dir+"/"
				+"SARS-CoV-2-Variants.txt";
		BufferedWriter bw = new BufferedWriter(new FileWriter( hapfil, false));
		String line ="";
		BufferedReader br = new BufferedReader(new FileReader(infil));
		while (( line = br.readLine()) != null) {
			if (!line.substring(0, 1).equals("#") ) {
				String[] tmp = line.split("\t");
				int index = Integer.parseInt( tmp[0].substring(1) );
				String freq = tmp[1];
				if (freq.length() > 7) {
					freq=freq.substring(0,6); 
				}
				if ( Double.parseDouble(freq) > this.Min_Hap_Freq ) {
					String hap= this.HaploArr.get(index);
					bw.write("SARS-CoV-2-Variant:\t"+ hap+"\t" + freq+"\n");
					String ss="SNP/INDEL";
//					System.out.println( this.Haplo2Genotype.get(hap));
					
					for (int i=0;i< this.Haplo2Genotype.get(hap).length();i++) {						
						if ( this.Haplo2Genotype.get(hap).substring(i,i+1).equals("1") ) {
							ss=ss+"\t" + this.VarArr.get(i);
						}
					}
					bw.write(ss+"\n");
				}
			}
		}
		br.close();
		bw.close();
            		
	}

	public void cmd(String command) throws IOException, InterruptedException {
		System.out.println(command);
		ProcessBuilder processBuilder = new ProcessBuilder();

		processBuilder.command("bash", "-c", command);

		try {
			Process process = processBuilder.start();

			InputStream inputStream = process.getInputStream();
			BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream));

			String line;
			while ((line = reader.readLine()) != null) {
				System.out.println(line);
			}

			process.waitFor();
		} catch (IOException | InterruptedException e) {
			e.printStackTrace();
		}
	}
}
