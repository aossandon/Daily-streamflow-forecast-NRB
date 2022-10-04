functions {
  
  matrix get_log_p_1cov_woint (vector covQ,real mu_mand, real sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand*covQ; 
     sigma = sigma_mand*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_1cov_woint_1_5 (vector covQ,real mu_mand, real sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand*(covQ .*sqrt(covQ)); 
     sigma = sigma_mand*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_1cov_woint_2 (vector covQ,real mu_mand, real sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand*(covQ .*covQ); 
     sigma = sigma_mand*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_1cov_woint_2_5 (vector covQ,real mu_mand, real sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand*(covQ .*covQ .*sqrt(covQ)); 
     sigma = sigma_mand*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_1cov_woint_3 (vector covQ,real mu_mand, real sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand*(covQ .*covQ .*covQ); 
     sigma = sigma_mand*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_1cov_woint_3_5 (vector covQ,real mu_mand, real sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand*(covQ .*covQ .*covQ .*sqrt(covQ)); 
     sigma = sigma_mand*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
  
   matrix get_log_p_1cov (vector covQ,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand[1] + mu_mand[2]*covQ; 
     sigma = sigma_mand[1] + sigma_mand[2]*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_1cov_1_5 (vector covQ,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand[1] + mu_mand[2]*(covQ .* sqrt(covQ)); 
     sigma = sigma_mand[1] + sigma_mand[2]*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_1cov_2 (vector covQ,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand[1] + mu_mand[2]*(covQ .* covQ); 
     sigma = sigma_mand[1] + sigma_mand[2]*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_1cov_2_5 (vector covQ,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand[1] + mu_mand[2]*(covQ .* covQ .* sqrt(covQ)); 
     sigma = sigma_mand[1] + sigma_mand[2]*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_1cov_3 (vector covQ,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand[1] + mu_mand[2]*(covQ .* covQ .* covQ); 
     sigma = sigma_mand[1] + sigma_mand[2]*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_1cov_3_5 (vector covQ,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand[1] + mu_mand[2]*(covQ .* covQ .* covQ .* sqrt(covQ)); 
     sigma = sigma_mand[1] + sigma_mand[2]*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_1cov_woint_sg (vector covQ,vector mu_mand, real sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand[1] + mu_mand[2]*covQ; 
     sigma = sigma_mand*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_1cov_woint_sg_1_5 (vector covQ,vector mu_mand, real sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand[1] + mu_mand[2]*(covQ .*sqrt(covQ));
     sigma = sigma_mand*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_1cov_woint_sg_2 (vector covQ,vector mu_mand, real sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand[1] + mu_mand[2]*(covQ .*covQ);
     sigma = sigma_mand*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_1cov_woint_sg_2_5 (vector covQ,vector mu_mand, real sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand[1] + mu_mand[2]*(covQ .*covQ .*sqrt(covQ));
     sigma = sigma_mand*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_1cov_woint_sg_3 (vector covQ,vector mu_mand, real sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand[1] + mu_mand[2]*(covQ .* covQ .*covQ);
     sigma = sigma_mand*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_1cov_woint_sg_3_5 (vector covQ,vector mu_mand, real sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand[1] + mu_mand[2]*(covQ .* covQ .*covQ .*sqrt(covQ));
     sigma = sigma_mand*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_1cov_woint_mu (vector covQ,real mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand*covQ; 
     sigma = sigma_mand[1] + sigma_mand[2]*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_1cov_woint_mu_1_5 (vector covQ,real mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand*(covQ .*sqrt(covQ)); 
     sigma = sigma_mand[1] + sigma_mand[2]*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_1cov_woint_mu_2 (vector covQ,real mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand*(covQ .*covQ); 
     sigma = sigma_mand[1] + sigma_mand[2]*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_1cov_woint_mu_2_5 (vector covQ,real mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand*(covQ .*covQ .*sqrt(covQ)); 
     sigma = sigma_mand[1] + sigma_mand[2]*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_1cov_woint_mu_3 (vector covQ,real mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand*(covQ .* covQ .*covQ); 
     sigma = sigma_mand[1] + sigma_mand[2]*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_1cov_woint_mu_3_5 (vector covQ,real mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand*(covQ .* covQ .*covQ .*sqrt(covQ)); 
     sigma = sigma_mand[1] + sigma_mand[2]*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_2cov (vector covQ, vector covP,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand[1] + mu_mand[2]*covQ + mu_mand[3]*covP; 
     sigma = sigma_mand[1] + sigma_mand[2]*covQ + sigma_mand[3]*covP;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_2cov_1_5 (vector covQ, vector covP,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand[1] + mu_mand[2]*(covQ .* sqrt(covQ)) + mu_mand[3]*covP; 
     sigma = sigma_mand[1] + sigma_mand[2]*covQ + sigma_mand[3]*covP;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_2cov_2 (vector covQ, vector covP,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand[1] + mu_mand[2]*(covQ .* covQ) + mu_mand[3]*covP; 
     sigma = sigma_mand[1] + sigma_mand[2]*covQ + sigma_mand[3]*covP;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_2cov_2_5 (vector covQ, vector covP,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand[1] + mu_mand[2]*(covQ .* covQ .* sqrt(covQ)) + mu_mand[3]*covP; 
     sigma = sigma_mand[1] + sigma_mand[2]*covQ + sigma_mand[3]*covP;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_2cov_3 (vector covQ, vector covP,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand[1] + mu_mand[2]*(covQ .* covQ .* covQ) + mu_mand[3]*covP; 
     sigma = sigma_mand[1] + sigma_mand[2]*covQ + sigma_mand[3]*covP;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_2cov_3_5 (vector covQ, vector covP,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand[1] + mu_mand[2]*(covQ .* covQ .* covQ .* sqrt(covQ)) + mu_mand[3]*covP; 
     sigma = sigma_mand[1] + sigma_mand[2]*covQ + sigma_mand[3]*covP;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
    
   
    matrix get_log_p_2cov_woi_sg_2cv (vector covQ, vector covP,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand[1] + mu_mand[2]*covQ + mu_mand[3]*covP; 
     sigma = sigma_mand[1]*covQ + sigma_mand[2]*covP;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_2cov_woi_sg_2cv_1_5 (vector covQ, vector covP,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand[1] + mu_mand[2]*(covQ .* sqrt(covQ)) + mu_mand[3]*covP; 
     sigma = sigma_mand[1]*covQ + sigma_mand[2]*covP;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_2cov_woi_sg_2cv_2 (vector covQ, vector covP,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand[1] + mu_mand[2]*(covQ .* covQ) + mu_mand[3]*covP; 
     sigma = sigma_mand[1]*covQ + sigma_mand[2]*covP;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_2cov_woi_sg_2cv_2_5 (vector covQ, vector covP,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand[1] + mu_mand[2]*(covQ .* covQ .* sqrt(covQ)) + mu_mand[3]*covP; 
     sigma = sigma_mand[1]*covQ + sigma_mand[2]*covP;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_2cov_woi_sg_2cv_3 (vector covQ, vector covP,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand[1] + mu_mand[2]*(covQ .* covQ .* covQ) + mu_mand[3]*covP; 
     sigma = sigma_mand[1]*covQ + sigma_mand[2]*covP;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_2cov_woi_sg_2cv_3_5 (vector covQ, vector covP,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand[1] + mu_mand[2]*(covQ .* covQ .* covQ .* sqrt(covQ)) + mu_mand[3]*covP; 
     sigma = sigma_mand[1]*covQ + sigma_mand[2]*covP;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_2cov_1sg (vector covQ, vector covP,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand[1] + mu_mand[2]*covQ + mu_mand[3]*covP; 
     sigma = sigma_mand[1] + sigma_mand[2]*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
    matrix get_log_p_2cov_1sg_1_5 (vector covQ, vector covP,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand[1] + mu_mand[2]*(covQ .* sqrt(covQ)) + mu_mand[3]*covP; 
     sigma = sigma_mand[1] + sigma_mand[2]*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
    matrix get_log_p_2cov_1sg_2 (vector covQ, vector covP,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand[1] + mu_mand[2]*(covQ .* covQ) + mu_mand[3]*covP; 
     sigma = sigma_mand[1] + sigma_mand[2]*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
  
   
   matrix get_log_p_2cov_1sg_2_5 (vector covQ, vector covP,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand[1] + mu_mand[2]*(covQ .* covQ .* sqrt(covQ)) + mu_mand[3]*covP; 
     sigma = sigma_mand[1] + sigma_mand[2]*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
    matrix get_log_p_2cov_1sg_3 (vector covQ, vector covP,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand[1] + mu_mand[2]*(covQ .* covQ .* covQ) + mu_mand[3]*covP; 
     sigma = sigma_mand[1] + sigma_mand[2]*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_2cov_1sg_3_5 (vector covQ, vector covP,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand[1] + mu_mand[2]*(covQ .* covQ .* covQ .* sqrt(covQ)) + mu_mand[3]*covP; 
     sigma = sigma_mand[1] + sigma_mand[2]*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_2cov_1sg_woi (vector covQ, vector covP,vector mu_mand, real sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand[1] + mu_mand[2]*covQ + mu_mand[3]*covP; 
     sigma = sigma_mand*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_2cov_1sg_woi_1_5 (vector covQ, vector covP,vector mu_mand, real sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand[1] + mu_mand[2]*(covQ .* sqrt(covQ)) + mu_mand[3]*covP; 
     sigma = sigma_mand*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
    matrix get_log_p_2cov_1sg_woi_2 (vector covQ, vector covP,vector mu_mand, real sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand[1] + mu_mand[2]*(covQ .* covQ) + mu_mand[3]*covP; 
     sigma = sigma_mand*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_2cov_1sg_woi_2_5 (vector covQ, vector covP,vector mu_mand, real sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand[1] + mu_mand[2]*(covQ .* covQ .* sqrt(covQ)) + mu_mand[3]*covP; 
     sigma = sigma_mand*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_2cov_1sg_woi_3 (vector covQ, vector covP,vector mu_mand, real sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand[1] + mu_mand[2]*(covQ .* covQ .* covQ) + mu_mand[3]*covP; 
     sigma = sigma_mand*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_2cov_1sg_woi_3_5 (vector covQ, vector covP,vector mu_mand, real sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu = mu_mand[1] + mu_mand[2]*(covQ .* covQ .* covQ .* sqrt(covQ)) + mu_mand[3]*covP; 
     sigma = sigma_mand*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_2cov_woi_b (vector covQ, vector covP,vector mu_mand, real sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu =  mu_mand[1]*covQ + mu_mand[2]*covP; 
     sigma = sigma_mand*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_2cov_woi_b_1_5 (vector covQ, vector covP,vector mu_mand, real sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu =  mu_mand[1]*(covQ .* sqrt(covQ)) + mu_mand[2]*covP; 
     sigma = sigma_mand*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_2cov_woi_b_2 (vector covQ, vector covP,vector mu_mand, real sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu =  mu_mand[1]*(covQ .* covQ) + mu_mand[2]*covP; 
     sigma = sigma_mand*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_2cov_woi_b_2_5 (vector covQ, vector covP,vector mu_mand, real sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu =  mu_mand[1]*(covQ .* covQ .* sqrt(covQ)) + mu_mand[2]*covP; 
     sigma = sigma_mand*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_2cov_woi_b_3 (vector covQ, vector covP,vector mu_mand, real sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu =  mu_mand[1]*(covQ .* covQ .* covQ) + mu_mand[2]*covP; 
     sigma = sigma_mand*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_2cov_woi_b_3_5 (vector covQ, vector covP,vector mu_mand, real sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu =  mu_mand[1]*(covQ .* covQ .* covQ .* sqrt(covQ)) + mu_mand[2]*covP; 
     sigma = sigma_mand*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
    matrix get_log_p_2cov_woi_musg (vector covQ, vector covP,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu =  mu_mand[1]*covQ + mu_mand[2]*covP; 
     sigma = sigma_mand[1]*covQ  + sigma_mand[2]*covP;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_2cov_woi_musg_1_5 (vector covQ, vector covP,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu =  mu_mand[1]*(covQ .* sqrt(covQ)) + mu_mand[2]*covP; 
     sigma = sigma_mand[1]*covQ + sigma_mand[2]*covP;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_2cov_woi_musg_2 (vector covQ, vector covP,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu =  mu_mand[1]*(covQ .* covQ) + mu_mand[2]*covP; 
     sigma = sigma_mand[1]*covQ + sigma_mand[2]*covP;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_2cov_woi_musg_2_5 (vector covQ, vector covP,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu =  mu_mand[1]*(covQ .* covQ .* sqrt(covQ)) + mu_mand[2]*covP; 
     sigma = sigma_mand[1]*covQ + sigma_mand[2]*covP;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_2cov_woi_musg_3 (vector covQ, vector covP,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu =  mu_mand[1]*(covQ .* covQ  .* covQ) + mu_mand[2]*covP; 
     sigma = sigma_mand[1]*covQ + sigma_mand[2]*covP;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_2cov_woi_musg_3_5 (vector covQ, vector covP,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu =  mu_mand[1]*(covQ .* covQ .* covQ .* sqrt(covQ)) + mu_mand[2]*covP; 
     sigma = sigma_mand[1]*covQ + sigma_mand[2]*covP;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_2cov_woi_mu (vector covQ, vector covP,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu =  mu_mand[1]*covQ + mu_mand[2]*covP; 
     sigma = sigma_mand[1] + sigma_mand[2]*covQ + sigma_mand[3]*covP;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_2cov_woi_mu_1_5 (vector covQ, vector covP,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu =  mu_mand[1]*(covQ .* sqrt(covQ)) + mu_mand[2]*covP; 
     sigma = sigma_mand[1] + sigma_mand[2]*covQ + sigma_mand[3]*covP;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_2cov_woi_mu_2 (vector covQ, vector covP,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu =  mu_mand[1]*(covQ .* covQ) + mu_mand[2]*covP; 
     sigma = sigma_mand[1] + sigma_mand[2]*covQ + sigma_mand[3]*covP;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_2cov_woi_mu_2_5 (vector covQ, vector covP,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu =  mu_mand[1]*(covQ .* covQ .* sqrt(covQ)) + mu_mand[2]*covP; 
     sigma = sigma_mand[1] + sigma_mand[2]*covQ + sigma_mand[3]*covP;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_2cov_woi_mu_3 (vector covQ, vector covP,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu =  mu_mand[1]*(covQ .* covQ .* covQ) + mu_mand[2]*covP; 
     sigma = sigma_mand[1] + sigma_mand[2]*covQ + sigma_mand[3]*covP;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_2cov_woi_mu_3_5 (vector covQ, vector covP,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu =  mu_mand[1]*(covQ .* covQ .* covQ .* sqrt(covQ)) + mu_mand[2]*covP; 
     sigma = sigma_mand[1] + sigma_mand[2]*covQ + sigma_mand[3]*covP;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_2cov_woi_mu_1csg (vector covQ, vector covP,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu =  mu_mand[1]*covQ + mu_mand[2]*covP; 
     sigma = sigma_mand[1] + sigma_mand[2]*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_2cov_woi_mu_1csg_1_5 (vector covQ, vector covP,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu =  mu_mand[1]*(covQ .* sqrt(covQ)) + mu_mand[2]*covP; 
     sigma = sigma_mand[1] + sigma_mand[2]*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_2cov_woi_mu_1csg_2 (vector covQ, vector covP,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu =  mu_mand[1]*(covQ .* covQ) + mu_mand[2]*covP; 
     sigma = sigma_mand[1] + sigma_mand[2]*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_2cov_woi_mu_1csg_2_5 (vector covQ, vector covP,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu =  mu_mand[1]*(covQ .* covQ .* sqrt(covQ)) + mu_mand[2]*covP; 
     sigma = sigma_mand[1] + sigma_mand[2]*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_2cov_woi_mu_1csg_3 (vector covQ, vector covP,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu =  mu_mand[1]*(covQ .* covQ .* covQ) + mu_mand[2]*covP; 
     sigma = sigma_mand[1] + sigma_mand[2]*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_2cov_woi_mu_1csg_3_5 (vector covQ, vector covP,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu =  mu_mand[1]*(covQ .* covQ .* covQ .* sqrt(covQ)) + mu_mand[2]*covP; 
     sigma = sigma_mand[1] + sigma_mand[2]*covQ;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   
   matrix get_log_p_2covmu_woi_Qcsg (vector covQ, vector covP,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu =  mu_mand[1] + mu_mand[2]*covQ + mu_mand[3]*covP; 
     sigma = sigma_mand[1] + sigma_mand[2]*covP;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_2covmu_woi_Qcsg_1_5 (vector covQ, vector covP,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu =  mu_mand[1] + mu_mand[2]*(covQ .* sqrt(covQ)) + mu_mand[3]*covP; 
     sigma = sigma_mand[1] + sigma_mand[2]*covP;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_2covmu_woi_Qcsg_2 (vector covQ, vector covP,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu =  mu_mand[1] + mu_mand[2]*(covQ .* covQ) + mu_mand[3]*covP; 
     sigma = sigma_mand[1] + sigma_mand[2]*covP;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
   matrix get_log_p_2covmu_woi_Qcsg_2_5 (vector covQ, vector covP,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu =  mu_mand[1] + mu_mand[2]*(covQ .* covQ .* sqrt(covQ )) + mu_mand[3]*covP; 
     sigma = sigma_mand[1] + sigma_mand[2]*covP;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
    matrix get_log_p_2covmu_woi_Qcsg_3 (vector covQ, vector covP,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu =  mu_mand[1] + mu_mand[2]*(covQ .* covQ .* covQ) + mu_mand[3]*covP; 
     sigma = sigma_mand[1] + sigma_mand[2]*covP;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
   
    matrix get_log_p_2covmu_woi_Qcsg_3_5 (vector covQ, vector covP,vector mu_mand, vector sigma_mand){
     matrix[rows(covQ) ,2] par;
     vector[rows(covQ)] mu;
     vector[rows(covQ)] sigma;
     int N;
     
     N=rows(covQ);
     mu =  mu_mand[1] + mu_mand[2]*(covQ .* covQ .* covQ .* sqrt(covQ )) + mu_mand[3]*covP; 
     sigma = sigma_mand[1] + sigma_mand[2]*covP;
     par[1:N,1] = (mu .*mu) ./ (sigma .* sigma);
     par[1:N,2] = mu ./ (sigma .* sigma);
     return par;
   }
} 

data {
  int<lower=0> N;
  int<lower=0> No;
  int<lower=0> S;
  int<lower=0> N_g1;
  int<lower=0> N_g2;
  int<lower=0> N_m1;
  int<lower=0> N_m2;
  int<lower=0> N_h1;
  int<lower=0> N_h2;
  int<lower=0> N_ho1;
  int<lower=0> N_ho2;
  int<lower=0> N_s1;
  int<lower=0> N_s2;
  
  int<lower=0> Np_g1;
  int<lower=0> Np_g2;
  int<lower=0> Np_m1;
  int<lower=0> Np_m2;
  int<lower=0> Np_h1;
  int<lower=0> Np_h2;
  int<lower=0> Np_ho1;
  int<lower=0> Np_ho2;
  int<lower=0> Np_s1;
  int<lower=0> Np_s2;
  matrix [No, S] y;
  matrix [No, S] covQ;
  matrix [No, S] covP;
  matrix [N, S] covQs;
  matrix [N, S] covPs;
  int indc[S];
  int ind_o[No];
  int ind_g1[N_g1];
  int ind_g2[N_g2];
  int ind_m1[N_m1];
  int ind_m2[N_m2];
  int ind_h1[N_h1];
  int ind_h2[N_h2];
  int ind_ho1[N_ho1];
  int ind_ho2[N_ho2];
  int ind_s1[N_s1];
  int ind_s2[N_s2];
  
  int indp_g1[Np_g1];
  int indp_g2[Np_g2];
  int indp_m1[Np_m1];
  int indp_m2[Np_m2];
  int indp_h1[Np_h1];
  int indp_h2[Np_h2];
  int indp_ho1[Np_ho1];
  int indp_ho2[Np_ho2];
  int indp_s1[Np_s1];
  int indp_s2[Np_s2];
  
  
}



parameters {
  vector<lower=0>[6] mu_garud;
  vector<lower=0>[6] sigma_garud;
  vector<lower=0>[6] mu_mand;
  vector<lower=0>[6] sigma_mand;
  vector<lower=0>[6] mu_han;
  vector<lower=0>[6] sigma_han;
  vector<lower=0>[6] mu_hos;
  vector<lower=0>[6] sigma_hos;
  vector<lower=0>[6] mu_san;
  vector<lower=0>[6] sigma_san;
  
 
}

model {
  matrix[No,2] tmp_gf;
  matrix[No,2] tmp_mf;
  matrix[No,2] tmp_hf;
  matrix[No,2] tmp_hof;
  matrix[No,2] tmp_sf;
  
  
  mu_garud ~ normal(0,10);
  sigma_garud ~ normal(0,10);
  mu_mand ~ normal(0,10);
  sigma_mand ~ normal(0,10);
  mu_han ~ normal(0,10);
  sigma_han ~ normal(0,10);
  mu_hos ~ normal(0,10);
  sigma_hos ~ normal(0,10);
  mu_san ~ normal(0,10);
  sigma_san ~ normal(0,10);
  
  // Garudeshwar
  tmp_gf[ind_g1,1:2] = get_log_p_2cov_2(covQ[ind_g1,indc[1]], covP[ind_g1,1],mu_garud[1:3], sigma_garud[1:3]);
  tmp_gf[ind_g2,1:2] = get_log_p_2cov_2(covQ[ind_g2,indc[1]], covP[ind_g2,1],mu_garud[4:6], sigma_garud[4:6]);
  
  // Mandleshwar
  tmp_mf[ind_m1,1:2] = get_log_p_2cov_2_5(covQ[ind_m1,indc[2]], covP[ind_m1,2],mu_mand[1:3], sigma_mand[1:3]);
  tmp_mf[ind_m2,1:2] = get_log_p_2cov_2_5(covQ[ind_m2,indc[2]], covP[ind_m2,2],mu_mand[4:6], sigma_mand[4:6]);
  
  //Handia
  tmp_hf[ind_h1,1:2] = get_log_p_2cov_1_5(covQ[ind_h1,indc[3]], covP[ind_h1,3],mu_han[1:3], sigma_han[1:3]);
  tmp_hf[ind_h2,1:2] = get_log_p_2cov_1_5(covQ[ind_h2,indc[3]], covP[ind_h2,3],mu_han[4:6], sigma_han[4:6]);
  
  //Hoshangabad
  tmp_hof[ind_ho1,1:2] = get_log_p_2cov_2(covQ[ind_ho1,indc[4]], covP[ind_ho1,4],mu_hos[1:3], sigma_hos[1:3]);
  tmp_hof[ind_ho2,1:2] = get_log_p_2cov_2(covQ[ind_ho2,indc[4]], covP[ind_ho2,4],mu_hos[4:6], sigma_hos[4:6]);
  
  //Sandiya
  tmp_sf[ind_s1,1:2] = get_log_p_2cov_2(covQ[ind_s1,indc[5]], covP[ind_s1,5],mu_san[1:3], sigma_san[1:3]);
  tmp_sf[ind_s2,1:2] = get_log_p_2cov_2(covQ[ind_s2,indc[5]], covP[ind_s2,5],mu_san[4:6], sigma_san[4:6]);
  
  
  
  target += gamma_lpdf(col(y,1) | col(tmp_gf,1), col(tmp_gf,2)) + gamma_lpdf(col(y,2) | col(tmp_mf,1), col(tmp_mf,2)) + gamma_lpdf(col(y,3) | col(tmp_hf,1), col(tmp_hf,2)) + gamma_lpdf(col(y,4) | col(tmp_hof,1), col(tmp_hof,2)) + gamma_lpdf(col(y,5) | col(tmp_sf,1), col(tmp_sf,2)); 
  
}


generated quantities{
  vector[No] logp_garud;
  vector[No] logp_man;
  vector[No] logp_han;
  vector[No] logp_hos;
  vector[No] logp_san;
  vector[No] log_lik;
  matrix<lower=0> [N, S] y_rep;
  matrix[N,2] tmp_gf;
  matrix[N,2] tmp_mf;
  matrix[N,2] tmp_hf;
  matrix[N,2] tmp_hof;
  matrix[N,2] tmp_sf;
  
  

  
  
  // for the year left out
  // Garudeshwar
  tmp_gf[indp_g1,1:2] = get_log_p_2cov_2(covQs[indp_g1,indc[1]], covPs[indp_g1,1],mu_garud[1:3], sigma_garud[1:3]);
  tmp_gf[indp_g2,1:2] = get_log_p_2cov_2(covQs[indp_g2,indc[1]], covPs[indp_g2,1],mu_garud[4:6], sigma_garud[4:6]);
  
  // Mandleshwar
  tmp_mf[indp_m1,1:2] = get_log_p_2cov_2_5(covQs[indp_m1,indc[2]], covPs[indp_m1,2],mu_mand[1:3], sigma_mand[1:3]);
  tmp_mf[indp_m2,1:2] = get_log_p_2cov_2_5(covQs[indp_m2,indc[2]], covPs[indp_m2,2],mu_mand[4:6], sigma_mand[4:6]);
  
  //Handia
  tmp_hf[indp_h1,1:2] = get_log_p_2cov_1_5(covQs[indp_h1,indc[3]], covPs[indp_h1,3],mu_han[1:3], sigma_han[1:3]);
  tmp_hf[indp_h2,1:2] = get_log_p_2cov_1_5(covQs[indp_h2,indc[3]], covPs[indp_h2,3],mu_han[4:6], sigma_han[4:6]);
  
  //Hoshangabad
  tmp_hof[indp_ho1,1:2] = get_log_p_2cov_2(covQs[indp_ho1,indc[4]], covPs[indp_ho1,4],mu_hos[1:3], sigma_hos[1:3]);
  tmp_hof[indp_ho2,1:2] = get_log_p_2cov_2(covQs[indp_ho2,indc[4]], covPs[indp_ho2,4],mu_hos[4:6], sigma_hos[4:6]);
  
  //Sandiya
  tmp_sf[indp_s1,1:2] = get_log_p_2cov_2(covQs[indp_s1,indc[5]], covPs[indp_s1,5],mu_san[1:3], sigma_san[1:3]);
  tmp_sf[indp_s2,1:2] = get_log_p_2cov_2(covQs[indp_s2,indc[5]], covPs[indp_s2,5],mu_san[4:6], sigma_san[4:6]);
  
  
  
  for (t in 1:N){
    
    if (t<=No){
      // Garudeshwar
      logp_garud[t] = gamma_lpdf(y[t,1] | tmp_gf[ind_o[t],1], tmp_gf[ind_o[t],2]);
      
      // Mandleshwar
      logp_man[t] = gamma_lpdf(y[t,2] | tmp_mf[ind_o[t],1], tmp_mf[ind_o[t],2]);
      
      //Handia
      logp_han[t] = gamma_lpdf(y[t,3] | tmp_hf[ind_o[t],1], tmp_hf[ind_o[t],2]);
    
      //Hoshangabad
      logp_hos[t] = gamma_lpdf(y[t,4] | tmp_hof[ind_o[t],1], tmp_hof[ind_o[t],2]);
    
      //Sandiya
      logp_san[t] = gamma_lpdf(y[t,4] | tmp_sf[ind_o[t],1], tmp_sf[ind_o[t],2]);
    }
    
    // Garudeshwar
    y_rep[t,1] = gamma_rng(tmp_gf[t,1], tmp_gf[t,2]);
      
    // Mandleshwar
    y_rep[t,2] = gamma_rng(tmp_mf[t,1], tmp_mf[t,2]);
    
    //Handia
    y_rep[t,3] = gamma_rng(tmp_hf[t,1], tmp_hf[t,2]);
    
    //Hoshangabad
    y_rep[t,4] = gamma_rng(tmp_hof[t,1], tmp_hof[t,2]);
    
    //Sandiya
    y_rep[t,5] = gamma_rng(tmp_sf[t,1], tmp_sf[t,2]);
    

  }
  

  log_lik = logp_garud + logp_man + logp_han + logp_hos + logp_san;
  
}

