function Autocorr(impurity)
if nargin==0
  impurity=[0.01,0.01,0.01,0.01]; %default impurities
end
[file,path]=uigetfile('*.csv','MultiSelect','on');
if iscell(file)
   fname=fullfile(path,file);     
else
   if file==0
      fprintf('user cancel selected\n');
      return
   else
      fname{1}=fullfile(path,file);
   end
end

%-----------------------------------------------------
tr={'C','N','D','O'}; %tracer symbles (do not change)
tracer={'13C','15N','2D','18O'};
% trA=1;  % specify 1st tracer A: 13C=1, 15N=2, 2D=3, 18O=4
% trB=3;  % specify 2nd tracer B
abundance=[0.0107,0.00364,0.00001,0.00187,0,0,0];

autogrouping=0;
warning('off','all')
%%--------start
for fnum=1:length(fname)
    T=readtable(fname{fnum},'readvariablename',true); 
    if isempty(T.Properties.VariableDescriptions)
        T.Properties.VariableDescriptions=T.Properties.VariableNames;
    end
    T=T(1:length(find([T.medMz]>0)),:); %cut empty rows.
    start_col=find(strcmp(T.Properties.VariableNames,'parent'))+1; %auto find start_rol
    sample_name=T.Properties.VariableNames(start_col:end)';
    grp_name=sample_name; %default grp_name, without grouping
    if autogrouping==1
         for i=1:length(sample_name)
             C=strsplit(sample_name{i},'_');
             if length(C)>1
                 grp_name{i,1}=sample_name{i}(1:length(sample_name{i})-length(C{end})-1);
             else
                 grp_name{i,1}=sample_name{i};
             end
         end 
    end
    grpHead=find(strcmp(T.isotopeLabel,'C12 PARENT')); 
    grpHead(end+1)=size(T,1)+1; %add a fake grphead(last)
    meta=[];
    for i=1:length(grpHead)-1        
            ids= grpHead(i):grpHead(i+1)-1; % rows ids for each metabolite
            T.metaGroupId(ids)=ones(1,length(ids))*i;
            A_sub=T(ids,:); %A_sub: data sheet for the selected metabolite ID
            dt = A_sub{:,start_col:end};
            try  %8/16/2020  added warning message for incorrect formula
                [~,~,tp]=formula2mass(A_sub.formula{1});
                
            catch
                fprintf(['check row#',num2str(ids(1)+1),' for errors in the formula name: ',A_sub.formula{1}],'Error detected!');
                return
            end
    
            lb=A_sub.isotopeLabel; 
            out=label_autodetect(lb); %auto detect C N H O
            if isempty(out)
                trA=0;
                trB=0;
                txt='no tracer';
                corr_abs=dt;
                corr_pct=ones(1,size(dt,2));
            elseif length(out)==1            
                trA=out(1);
                trB=out(1);
                txt=tracer{trA};
                A_num=tp(trA);  % A num
                ab_A=abundance(trA);
                impurity_A=impurity(trA);
                counts = AB_getcounts(lb,tr{trA},tr{trA});
                [fulldt,reduced_idx] = AB_getfulldt(dt,counts,A_num,0);
                [corr_abs,corr_pct]=isocorr_A(fulldt,A_num,ab_A,impurity_A);
                corr_abs=corr_abs(reduced_idx,:); %shorttable for output
                corr_pct=corr_pct(reduced_idx,:); %shorttable for output            
           
            elseif length(out)==2            
                trA=out(1);
                trB=out(2);
                txt=[tracer{trA},'&',tracer{trB}];
                A_num=tp(trA);  % A num
                B_num=tp(trB);  % B num
                ab_A=abundance(trA);
                ab_B=abundance(trB); 
                impurity_A=impurity(trA);
                impurity_B=impurity(trB);
                counts = AB_getcounts(lb,tr{trA},tr{trB});            
                [fulldt,reduced_idx] = AB_getfulldt(dt,counts,A_num,B_num);
                [corr_abs,corr_pct]=isocorr_AB(fulldt,A_num,B_num,ab_A,ab_B,impurity_A,impurity_B);
                corr_abs=corr_abs(reduced_idx,:); %shorttable for output
                corr_pct=corr_pct(reduced_idx,:); %shorttable for output            
            else
                fprintf("error: more than 2 labeled terms");
            end
            %store individuals 
            meta(i).ID=i;      
            meta(i).name=A_sub.compound{1};
            meta(i).formula=A_sub.formula{1};
            meta(i).mz=A_sub.medMz(1);
            meta(i).rt=A_sub.medRt(1);
            meta(i).corr_abs=corr_abs; %shorttable for output
            meta(i).corr_pct=corr_pct; %shorttable for output
            meta(i).corr_tic=sum(corr_abs,1);
    
            fprintf([num2str(i),'/',num2str(length(grpHead)-1),': '])           
            fprintf(txt);
            fprintf([' -- ',meta(i).name])
            [max_ppm,idx]=max(abs([A_sub.ppmDiff])); 
            errline=find([A_sub.ppmDiff]==max_ppm)+ids(1);
            if max_ppm>5
                fprintf(2,[' -- warning: line',num2str(errline),', ppm error larger than expected:',num2str(max_ppm,'%.2f'),'\n']);
            else
                fprintf(' -- successful\n')
            end
           
    end 
    % -------- concat
    cat_abs=[];cat_pct=[];
    for i=1:length(meta)  
     cat_abs=[cat_abs;meta(i).corr_abs];  %concatenate for csv output
     cat_pct=[cat_pct;meta(i).corr_pct];  
    end
    
    % --------- put together into tables
    A_part1=T(:,1:start_col-1);
    A_part2=T(:,start_col:end);
    
    A_part2{:,:}=cat_abs;
    A_corr_abs=[A_part1,A_part2];
    
    A_part2{:,:}=cat_pct;
    A_corr_pct=[A_part1,A_part2];
    
    % make 3rd table as total ion/////////////
    A_corr_total=[];
    for i=1:length(meta)
      A_corr_total{i,1}=meta(i).ID;
      A_corr_total{i,2}=meta(i).name;
      A_corr_total{i,3}=meta(i).formula;
      for j=1:length(meta(i).corr_tic)
          A_corr_total{i,3+j}=meta(i).corr_tic(j);
      end
    end
    A_corr_total=cell2table(A_corr_total);
    %A_corr_total.Properties.VariableNames=[{'ID','Name','formula'},T.Properties.VariableNames(start_col:end)];
    A_corr_total.Properties.VariableNames=[{'ID','Name','formula'},T.Properties.VariableDescriptions(start_col:end)];
    
   % end 3rd table  ////////////////////////
    
    %save 
    
    [filepath,name,~] = fileparts(fname{fnum});

    T.Properties.VariableNames=T.Properties.VariableDescriptions;
    A_corr_pct.Properties.VariableNames=T.Properties.VariableDescriptions;
    A_corr_abs.Properties.VariableNames=T.Properties.VariableDescriptions;
    
    fname_all=fullfile(filepath,[name,'_cor','.xlsx']);
    writetable(T,fname_all,'Sheet','original');
    writetable(A_corr_pct,fname_all,'Sheet','cor_pct');
    writetable(A_corr_abs,fname_all,'Sheet','cor_abs');
    writetable(A_corr_total,fname_all,'Sheet','total');
    fprintf(['file #',num2str(fnum),': ',name,'------ done!\n'])
end




