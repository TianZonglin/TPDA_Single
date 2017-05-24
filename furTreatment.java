import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Vector;
/*
 * This program is the re-processing of the Apriori algorithm's results.
 * The main content is the consolidation of association rules.
 */
public class furTreatment {
	

	@SuppressWarnings({ "rawtypes", "unchecked" })
	public static void main(String args[]) throws IOException{
		
 
		 int r = Integer.valueOf(args[0])+1; 		//the number of rules after Aprori
		 String f = args[1];						//the path of input file
		 String fout = args[2];						//the path of output file
		 String di = ","; 
		 
		 File file = new File(fout);
		 //If file doesn't exists,create it.
		  if (!file.exists()) {
			  file.createNewFile();
		  }
		  else{
			  try {
				  FileWriter fw = new FileWriter(fout,false);
				  BufferedWriter bw = new BufferedWriter(fw);
				  
		          bw.write("");
		          bw.close();
		      } catch (Exception e) {
		          e.printStackTrace();
		      }
		  }
		
		try{					
			BufferedReader b = new BufferedReader(new FileReader(f));
			FileWriter fw = new FileWriter(fout,true);
			BufferedWriter bw = new BufferedWriter(fw);
			String line;
			int j = 0;
			Vector Data[] = new Vector[r-1];
			System.out.println("Start");
			b.readLine(); 
			
			//Read all data from the file
			while((line = b.readLine()) != null){
				line = line.replace("\"", ""); 
				line = line.replace("{", ""); 
				line = line.replace("}", ""); 
				line = line.replace(" ", ""); 
				// All rules use "=>" on both sides to divide into two parts
				String rules[] = line.split("=>"); 
				String front[] = rules[0].split(di);
				
				//Remove the 'many-to-one' rules	
              	if(!(front[1].length() > 20)){
              		if(Data[j] == null){
              			Data[j] = new Vector<String>();
              		}
              		
              		String head[] = front[1].split("[*]");//Split again,split the number.
              		
              		Data[j].addAll(Arrays.asList(head[0]));//0 gene name
              		Data[j].addAll(Arrays.asList(head[1]));//1 gene part number
              		
              		//Do with laster part
              		String item[] = rules[1].split(","); 
              		//Remove the 'one-to-many' rules
              		if(!(item.length > 4)){
              			
              			String tail[] = item[0].split("[*]");//Split again,split the number.
                  		
                  		Data[j].addAll(Arrays.asList(tail[0]));//2 laster gene name
                  		Data[j].addAll(Arrays.asList(tail[1]));//3 laster gene part number
                  		
                      	Data[j].addAll(Arrays.asList(item[1]));//4 support value
                      	j++;
              		}
              		else{
              			Data[j].clear();
      					j--; 
              		}
              	}
			} 
			
			//Change the form 'B->A' to 'A->B'
			for(int i = 0;i < j-1;i++){ 
	        	  //System.out.println(i);
				  String t;				  
				  if((Data[i].get(0).toString()).compareTo(Data[i].get(2).toString()) > 0){
					  t = Data[i].get(0).toString();//gene
					  Data[i].set(0,Data[i].get(2).toString());
					  Data[i].set(2,t);
					  
					  t = Data[i].get(1).toString();//number
					  Data[i].set(1,Data[i].get(3).toString());
					  Data[i].set(3,t);
				  }
	      	}
			
			//Delete duplicate
			for(int k = 0;k < j-1;k++){
				for(int i = k+1;i < j;i++){
		        	  //System.out.println(i);
					  if((Data[k].get(0).toString()).equals(Data[i].get(0).toString()) 
						  && (Data[k].get(2).toString()).equals(Data[i].get(2).toString())
							&&(Data[k].get(1).toString()).equals(Data[i].get(1).toString()) 
								&& (Data[k].get(3).toString()).equals(Data[i].get(3).toString())){
						  Data[k].set(0,"$");
						  Data[k].set(2,"$");
					  }
		      	}
			}

			//The first step of sort : insert the sort by the front.
			Vector<String> t1 = new Vector<String>();
			int temp1;
	        for(int k1 = 0 ;k1 <= j-1;k1++){ 
	        	temp1 = k1;
	      		for(int k2 = k1;k2 <= j-1;k2++){
	      			if(((Data[k2].get(0)).toString()).compareTo((Data[temp1].get(0)).toString()) > 0){
      					temp1 = k2;
      				}	
	      		}
	      		if(temp1 != k1){
	      			t1 = Data[k1];
	      			Data[k1] = Data[temp1];
	      			Data[temp1] = t1;
	      		}
	      	 }
	         //The first step of sort : When the front is equal, insert the sort by the rear.
	         int old = 0;
			 for(int k3 = 0;k3 <= j-1;k3++){
				 if(!(Data[old].get(0).toString()).equals(Data[k3].get(0).toString()) || (k3 == j-1)){//Make sure the end is right.
					 Vector<String> t2 = new Vector<String>();
					 int temp2;
				     for(int k1 = old ;k1 < k3;k1++){
				        temp2 = k1;
				      	for(int k2 = k1;k2 < k3;k2++){
				      		if(((Data[k2].get(2)).toString()).compareTo((Data[temp2].get(2)).toString()) > 0){
				      			temp2 = k2;
				      		}
				      	}
				      	if(temp2 != k1){
				      		t2 = Data[k1];
				      		Data[k1] = Data[temp2];
				      		Data[temp2] = t2;
				      	}
				      }
				      old = k3;
				 } 			 
			 }

			 //merge data
			 int anchor = 0;
			 for(int k2 = 0;k2 <= j-1;k2++){
				 if(!(Data[anchor].get(0).toString()+Data[anchor].get(2).toString()).equals(Data[k2].get(0).toString()+Data[k2].get(2).toString()) || (k2 == j-1)){
						double refer = 0;
						for(int k1 = anchor ;k1 < k2;k1++){							 
							 refer = refer + Double.parseDouble(Data[k1].get(4).toString());
							 if(k1 > anchor){
								Data[k1].clear();
								Data[k1].addAll(Arrays.asList("$"));
							}
						}
						Data[anchor].set(4,refer);
						Data[anchor].set(1,"sum");
						Data[anchor].set(3,"sum");
						anchor = k2;
					 }	 			 
			 }
          	
			 //output results
	         bw.write("h_rule,h_no,t_rule,t_no,reference" + "\r\n");
	         for(int i = 0;i < j;i++){ 
	        	      		
	        	 if(!(Data[i].get(0).toString()).equals("$")){
	        	     bw.write(Data[i].toString().replace("[", "").replace("]", "").replace(" ", "") + "\r\n");
	        	 }
	         }
          	
			System.out.println("End");
			System.out.println("Next SortbyRules.r");
			bw.close();
			b.close();
		}catch (Exception e) {
	        e.printStackTrace();
		}
	}
}
