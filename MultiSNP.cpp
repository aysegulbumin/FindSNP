//
// Created by aysegull on 10/5/18.
//


#include <sstream>
#include <map>
#include <algorithm>
#include <set>
#include <functional>
#include <cmath>
#include <iostream>
#include <fstream>
#include <chrono>

#include "MultiSNP.h"

using namespace std;
using namespace std::chrono;

map< string, string > Dictionary_;
map <string, map<string,int>> from_where;
int global_count=0;
int no_change_in_sequence_count=0;
int change_in_sequence=0;
int flag_changed=0;
int flag_changed_count=0;
int flag_notchanged_count=0;
map <string, int> gene_count;
map <string, int> gene_count2;
int global_total_mut=0;
int global_btw_ss=0;

class time_point;

void Initialize2()
{

    Dictionary_["TTT"]="phe";
    Dictionary_["TTC"]="phe";
    Dictionary_["TTA"]="leu";
    Dictionary_["TTG"]="leu";
    Dictionary_["CTT"]="leu";
    Dictionary_["CTC"]="leu";
    Dictionary_["CTA"]="leu";
    Dictionary_["CTG"]="leu";
    Dictionary_["ATT"]="ile";
    Dictionary_["ATC"]="ile";
    Dictionary_["ATA"]="ile";
    Dictionary_["ATG"]="met";
    Dictionary_["GTT"]="val";
    Dictionary_["GTC"]="val";
    Dictionary_["GTA"]="val";
    Dictionary_["GTG"]="val";
    Dictionary_["TCT"]="ser";
    Dictionary_["TCC"]="ser";
    Dictionary_["TCA"]="ser";
    Dictionary_["TCG"]="ser";
    Dictionary_["CCT"]="pro";
    Dictionary_["CCC"]="pro";
    Dictionary_["CCA"]="pro";
    Dictionary_["CCG"]="pro";
    Dictionary_["ACT"]="thr";
    Dictionary_["ACC"]="thr";
    Dictionary_["ACA"]="thr";
    Dictionary_["ACG"]="thr";
    Dictionary_["GCT"]="ala";
    Dictionary_["GCC"]="ala";
    Dictionary_["GCA"]="ala";
    Dictionary_["GCG"]="ala";
    Dictionary_["TAT"]="tyr";
    Dictionary_["TAC"]="tyr";
    Dictionary_["CAT"]="his";
    Dictionary_["CAC"]="his";
    Dictionary_["CAG"]="gln";
    Dictionary_["CAA"]="gln";
    Dictionary_["AAT"]="asn";
    Dictionary_["AAC"]="asn";
    Dictionary_["AAA"]="lys";
    Dictionary_["AAG"]="lys";
    Dictionary_["GAT"]="asp";
    Dictionary_["GAC"]="asp";
    Dictionary_["GAA"]="glu";
    Dictionary_["GAG"]="glu";
    Dictionary_["TGT"]="cys";
    Dictionary_["TGC"]="cys";
    Dictionary_["TGG"]="trp";
    Dictionary_["CGT"]="arg";
    Dictionary_["CGC"]="arg";
    Dictionary_["CGA"]="arg";
    Dictionary_["CGG"]="arg";
    Dictionary_["AGT"]="ser";
    Dictionary_["AGC"]="ser";
    Dictionary_["AGA"]="arg";
    Dictionary_["AGG"]="arg";
    Dictionary_["GGT"]="gly";
    Dictionary_["GGC"]="gly";
    Dictionary_["GGA"]="gly";
    Dictionary_["GGG"]="gly";
    Dictionary_["TGA"]="STOP";
    Dictionary_["TAA"]="STOP";
    Dictionary_["TAG"]="STOP";
}

vector<string> convertbytabb(string input)
{
    vector<char> v;

    istringstream ss(input);
    string token;
    vector<string> line_vector;
    while(std::getline(ss, token, '\t'))
    {
        line_vector.push_back(token);
    }
    return line_vector;
}

void findAllOccurances_(std::vector<size_t> & vec, std::string data, std::string toSearch)
{
    // Get the first occurrence
    size_t pos = data.find(toSearch);

    // Repeat till end is reached
    while( pos != std::string::npos)
    {
        // Add position to the vector
        vec.push_back(pos);

        // Get the next occurrence from the current position
        pos =data.find(toSearch, pos + toSearch.size());
    }
}

void AminoAcidSequence2 (string original, string mutated,string genename) {
 //This part is implemented in case the starting and stopping codons do not matter - so not used for now
    string sequence="";
    string mutated_sequence="";
    for(int i=0;i<original.size()-original.size()%3;i=i+3) {
    //No matter whether starting codon and stopping codon, always build the entire sequence
        string one(1, original[i]);
        string two(1, original[i+1]);
        string three(1, original[i+2]); // The codon from the original -gene string is this one

        string one_mutated(1,mutated[i]);
        string two_mutated(1, mutated[i+1]);
        string three_mutated(1, mutated[i+2]); // The codon from the mutated string is this one.

        string codon="";
        string codon_mutated="";


        codon=one+two+three;// The codon from the original -gene string is this one

        codon_mutated=one_mutated+two_mutated+three_mutated;// The codon from the mutated string is this one.
        sequence=sequence+Dictionary_[codon]+" ";
        mutated_sequence=mutated_sequence+Dictionary_[codon_mutated]+" ";

    }
    if(sequence.compare(mutated_sequence)!=0)
    {
        if ( gene_count2.find(genename) == gene_count2.end() ) {
            // not found
            gene_count2[genename]=1;
        }
        else
        {
            gene_count2[genename]++;
        }
         //They are different
        flag_changed_count++;
        cout <<genename<<endl;
        cout<<sequence<<endl;
        cout << mutated_sequence<<endl;
    }
    else
    {
        flag_notchanged_count++;
    }

}
int AminoAcidSequence (string original, string mutated,string genename, int total_mut, ofstream& output_amino) {

    vector <int> indices;
    vector<vector <int> > allindices;
    vector <string> normal_ones;
    vector<string> mut_ones;

    string sequence="";
    string mutated_sequence="";
    size_t found = original.find("ATG"); //The first position that you found the ATG which is the starting codon.
    int j=static_cast<int>(found);

    for(int i=j+3;i<original.size();i=i+3) {
        //START FROM WHERE THE ATG starting codon ends.
        //and read 3 at a time and iterate by adding 3 each time.

        string one(1, original[i]);
        string two(1, original[i+1]);
        string three(1, original[i+2]); // The codon from the original -gene string is this one

        string one_mutated(1,mutated[i]);
        string two_mutated(1, mutated[i+1]);
        string three_mutated(1, mutated[i+2]); // The codon from the mutated string is this one.

        string codon="";
        string codon_mutated="";


        codon=one+two+three;// The codon from the original -gene string is this one

        codon_mutated=one_mutated+two_mutated+three_mutated;// The codon from the mutated string is this one.


        if(codon.compare("TAA")==0 || codon.compare("TAG")==0 || codon.compare("TGA")==0 )
        {
            //Stop codon found
            //Next index is found


            if(sequence.compare(mutated_sequence)!=0)
            {
                //cout<<"Normal:"<<sequence<<endl;
                //cout<<"Mutated:"<<mutated_sequence<<endl;
                //If they are different
                normal_ones.push_back(sequence);
                mut_ones.push_back(mutated_sequence);
                allindices.push_back(indices);
                change_in_sequence++;
                flag_changed=1;

            }
            else { no_change_in_sequence_count++;}

            std::string str2 = original.substr (i+3,original.size());
            std::size_t found = str2.find("ATG");
            i=static_cast<int>(found)+i+3;
            sequence="";
            mutated_sequence="";
            indices.clear();
        }
        else
        {
            sequence=sequence+Dictionary_[codon]+" ";
            mutated_sequence=mutated_sequence+Dictionary_[codon_mutated]+" ";
            if(Dictionary_[codon].compare(Dictionary_[codon_mutated])!=0)
            {
                indices.push_back(i);
                if(from_where.find(Dictionary_[codon])==from_where.end())
                {
                    map<string,int> temp_;
                    temp_[Dictionary_[codon_mutated]]=1;
                    from_where[Dictionary_[codon]]=temp_;
                }
                else
                {
                    if(from_where[Dictionary_[codon]].find(Dictionary_[codon_mutated])==from_where[Dictionary_[codon]].end())
                    {

                        from_where[Dictionary_[codon]][Dictionary_[codon_mutated]]=1;
                    }
                    else
                    {
                        from_where[Dictionary_[codon]][Dictionary_[codon_mutated]]+=1;
                    }
                }
            }

        }
    }
    if(flag_changed==1)
    {
        flag_changed_count++;
        flag_changed=0;
        if ( gene_count2.find(genename) == gene_count2.end() ) {
            // not found
            gene_count2[genename]=1;
        }
        else
        {
            gene_count2[genename]++;
        }
    }
    else if(flag_changed==0)
    {
        flag_notchanged_count++;
    }
    //cout<<  "Total Mut"<<total_mut<<endl;
    //cout<< "Between Start and Stop"<<total_mut-mut_ones.size()<<endl;
    global_total_mut+=total_mut;
    global_btw_ss+=total_mut-mut_ones.size();
    for(int l=0;l<mut_ones.size();l++)
    {

        output_amino<<genename<<endl;
        //cout<<"Size:"<<allindices[l].size()<<endl;
        for(auto k:allindices[l])
        {
            output_amino<<k<<" ";
        }
        output_amino<<endl;
        global_count++;
        //cout<<genes[l]<<endl;
        output_amino<<normal_ones[l]<<endl; //Actual aminoacid sequence.
        output_amino<<mut_ones[l]<<endl;    //Aminoacid sequence after the mutation.
        output_amino<<endl;
    }


    return mut_ones.size();

}

int MultiSNP(int argc, char** argv)
{
    //for CPU time
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    Initialize2();
    if(argc<7)
    {
        cout<< "Please give all the parameters.\nRespectively: \n \t Tab Delimited SAM file, Reference fasta file, threshold , SNP details output file, Amino acid sequence Details output file , Gene Details output file"<<endl;
        return 0;
    }
    int threshold=stoi(argv[3]);
    //AminoAcidSequence("ATGCAGTAAAAAATGAAGAAGTAATACATGATATAA","ATGAAATAAAAAATGAACAAGTAATAAATGAAATAA","GeneName");
    ifstream myfile;
    string line;
    myfile.open(argv[1]);
    ofstream SNPDetails;
    SNPDetails.open(argv[4]);
    ofstream AminoAcidSeqDetails;
    AminoAcidSeqDetails.open(argv[5]);
    map <string, int > flags;
    map <pair<string,string>, vector<vector<string> > > gene_Specific;
    pair <string, pair<string,string> > read_ref_flag  ;
    map <pair <string,string>, pair<map<string,int>,vector<pair<string,int>>>> gene; // key is gene name and position and the values were the list of mutations and the number of occurences of that mutation
    //Now the value is pair of string vector and vector od string and int pairs .
    map <pair <string,string>, string> gene2;
    read_ref_flag.first="";
    read_ref_flag.second.first="";
    read_ref_flag.second.second="";
    vector<string> bases={"A","C","G","T"};
    map< pair<string,string> , pair<vector<int>, vector<string>> > gene_read;
    //gene name, read name  as a key and for each read and gene pair we have the position of the mutation and the particular mutation that occured in that position
    map< string, map< pair<vector<int>, vector<string>>, int>> multisnp_count;

    if (myfile.is_open()) {

        while (getline(myfile, line)) {
            if (line[0] != '#') {

                vector<string> myvec = convertbytabb(line);
                if (stoi(myvec[1]) == 117 || stoi(myvec[1]) == 181 || stoi(myvec[1]) == 69 || stoi(myvec[1]) == 133 ||stoi(myvec[1]) == 141 || stoi(myvec[1]) == 77) {
                    //skip
                    //These are the flags for the unmapped reads

                } else {
                    if (myvec[4] != myvec[7]) {

                        //INDICATES MUTATION

                        if ((std::find(std::begin(bases), std::end(bases), myvec[4]) != std::end(bases)) &&
                            (std::find(std::begin(bases), std::end(bases), myvec[7]) != std::end(bases))) {
                            //INDICATES NO DELETION NO INSERTION

                            pair <string,string> temp_gene_read;
                            temp_gene_read.first=myvec[2]; //gene name
                            temp_gene_read.second=myvec[0]; //read name
                            if ( gene_read.find(temp_gene_read) == gene_read.end() ) {
                                // not found
                                vector<int> position_of_gene;
                                vector<string> mutation;
                                pair<vector<int>, vector<string>> PositionAndMutation;
                                position_of_gene.push_back(stoi(myvec[6]));
                                mutation.push_back(myvec[7]+"_"+myvec[4]);
                                PositionAndMutation.first=position_of_gene;
                                PositionAndMutation.second=mutation;
                                gene_read[temp_gene_read];
                                gene_read[temp_gene_read]=PositionAndMutation;


                            } else {
                                // found
                                if (std::find(gene_read[temp_gene_read].first.begin(),gene_read[temp_gene_read].first.end(),stoi(myvec[6]))!=gene_read[temp_gene_read].first.end())
                                {
                                    //that particular position and mutation pair is found already in the structure
                                }
                                else
                                {
                                    gene_read[temp_gene_read].first.push_back(stoi(myvec[6]));
                                    gene_read[temp_gene_read].second.push_back(myvec[7]+"_"+myvec[4]);

                                }

                            }
                        }


                    }

                }
            }

        }
        myfile.close();
    }

    map <int,int> counts;

    for(auto elem:gene_read)
    {
        if ( multisnp_count.find(elem.first.first) == multisnp_count.end() ) {
            // not found
            map< pair<vector<int>, vector<string>>, int> temp_map;
            temp_map[elem.second];
            temp_map[elem.second]=1;

            //multisnp_count[elem.first.first];
            multisnp_count[elem.first.first]=temp_map;


        } else {
            // found
            multisnp_count[elem.first.first][elem.second]++;

        }

    }


    ifstream myfile2;
    myfile2.open(argv[2]);
    string line2;
    string genename;
    string genesequence;
    map<string,string> name_sequence;
    if (myfile2.is_open()) {

        while (getline(myfile2, line2)) {
            if(line2[0]=='>')
            {
                genename=line2.substr(1, line2.size() - 1);
            }
            else
            {
                genesequence=line2;
                name_sequence[genename];
                name_sequence[genename]=genesequence;
              //  cout<<name_sequence[genename]<<endl;
            }
        }
        myfile2.close();
    }

    int k=0;
    cout<<  "SNP Details and Aminoacid Sequence Changes are being  calculated..."<<endl;
    for(auto e:multisnp_count)
    {
        for(auto e2:e.second)
        {
            if(e2.second>=threshold) //This mutation list and position list pair is seen in more than 10 reads if threshold is 10.
            {

                if(e2.first.first.size()>=1) //If e2.first.first.size()>=1 there is more than one element in the list of mutations, if it is equal to 1 this means there is only one SNP and it already occurs in more than 10 reads in the same position.
                {
                    //MORE THAN 2 SNP //2 is a parameter

                    SNPDetails<<e2.second<<" Reads have the following mutation in the gene "<< e.first<<endl;
                   // cout<<name_sequence[e.first]<<endl;

                    if ( gene_count.find(e.first) == gene_count.end() ) {
                        // not found
                        gene_count[e.first]=1;
                    }
                    else
                    {
                        gene_count[e.first]++;
                    }


                    string mutated=name_sequence[e.first];
                    for(int i=0;i<e2.first.first.size(); i++)
                    {
                        mutated[e2.first.first[i]-1]=e2.first.second[i][2]; //create the mutated sequence from the mutations that we have. The mutation was saved in form of From_To. So thats why the second index is considered.
                    }
                    //AminoAcidSequence("ATGCAGTAAAAAATGAAGAAGTAATACATGATATAA","ATGAAATAAAAAATGAACAAGTAATAAATGAAATAA","GeneName");


                    SNPDetails<<e2.first.first.size()<<" Mutations per read" <<endl;
                    SNPDetails<<"\tPosition in the Gene and The Mutation"<<endl;
                    k++;
                    for(int i=0;i<e2.first.first.size(); i++)
                    {
                       // cout<<name_sequence[e.first][e2.first.first[i]-1]<<endl;
                        SNPDetails<<"\t"<<e2.first.first[i]<<"\t"<<e2.first.second[i]<<endl;
                    }

                   int num_of_mut_btw_startstop=AminoAcidSequence(name_sequence[e.first],mutated,e.first,e2.first.first.size(),AminoAcidSeqDetails);
                   if(e2.first.first.size()- num_of_mut_btw_startstop>0)
                   {
                      cout<< "All"<< e2.first.first.size()<<endl;
                      cout<< "BTW start  and stop"<< num_of_mut_btw_startstop<<endl;
                     cout<< "Not btw start  and stop"<< e2.first.first.size()- num_of_mut_btw_startstop<<endl;
                     cout<<endl;
                   }
                  // AminoAcidSequence2(name_sequence[e.first],mutated,e.first);
                }


            }
        }
    }
    AminoAcidSeqDetails.close();
    SNPDetails.close();
    SNPDetails.close();
   // cout<<"In total "<<k<<endl;


    //int cou=0;
   /* for(auto t:from_where)
    {
        cout<<"From "<<t.first<<" ";
        for(auto l:from_where[t.first])
        {
            cout<<"\tTo ->"<<l.first<<" "<<l.second<<" Times"<<endl;
            cou=cou+l.second;
        }
    }*/
    //cout<<"Total changes"<<endl;

    //cout<<cou<<endl;

    /*cout<<  "Global count: "<<global_count<<endl;
    cout<<  "Change in sequence in this many sequences"<<change_in_sequence<<endl;
    cout<<  "No Change in sequence in this many sequences"<<no_change_in_sequence_count<<endl;
    */
    cout<<  "Gene Details file is being created..."<<endl;
    ofstream GeneDetails;
    GeneDetails.open(argv[6]);
    GeneDetails<<"There are "<<gene_count.size()<< " genes that have mutations."<<endl;
    int sum=0;
    GeneDetails<<"Gene\tCount"<<endl;
    for(auto elem: gene_count)
    {
        sum=sum+elem.second;
        GeneDetails<< elem.first<<"\t"<<elem.second<<endl;
    }
    GeneDetails<<"Adding up to "<<sum <<endl;

    GeneDetails<<"\n\n\n"<<endl;


    GeneDetails<<"There are "<<gene_count2.size()<< " genes that have change visible in the aminoacid sequence."<<endl;
    GeneDetails<<"Gene\tCount"<<endl;
    sum=0;
    for(auto elem: gene_count2)
    {
        sum=sum+elem.second;
        GeneDetails<<elem.first<<"\t"<<elem.second<<endl;
    }
    GeneDetails<<"Adding up to "<<sum <<endl;

    GeneDetails.close();
    //THIS PART WAS TO SEE THE OUTPUT OF THE FIRST STEP , GENE NAME, READ NAME AND ALL POSITIONS AND MUTATIONS.
    /*
    for(auto elem: gene_read)
    {
        if ( counts.find(elem.second.first.size()) == counts.end() ) {
            //not found
            counts[elem.second.first.size()];
            counts[elem.second.first.size()]=1;

        } else {
            // found
            counts[elem.second.first.size()]+=1;

        }

        cout<<  elem.second.first.size()<<endl;
        cout<<  elem.second.second.size()<<endl;

        cout<<"Gene Name:"<<elem.first.first<<"\tRead Name:"<<elem.first.second<<endl;

        cout<<"\t"<<"Mutation"<<"\t"<<"Position"<<endl;

        for(int i=0;i<elem.second.first.size(); i++)
        {
            cout<<"\t"<<elem.second.second[i]<<"\t"<<elem.second.first[i]<<endl;
        }

    }*/
    //THIS PART IS TO SEE HOW MANY ONLY K MUTATED READS EXIST K  CHANGES BETWEEN 1 AND 179 .
   /* cout<<"************************************************************************"<<endl;

    for(auto elem: counts)
    {
        cout<<elem.first<<"\t-->"<<elem.second<<endl;
    }
    */
    /*cout << "Actual change"<<endl;
    cout<<flag_changed_count<<endl;
    cout<<  "Actual no change"<<endl;
    cout<<flag_notchanged_count<<endl;

    cout<<"Global total mut"<<endl;
    cout<<global_total_mut<<endl;
    cout<<"Between start and stop"<<endl;
    cout<<global_btw_ss<<endl;
    */

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<seconds>( t2 - t1 ).count();
    cout<<"CPU Time: "<<duration<<" seconds. \n";
    cout<<"Done!"<<endl;
        return 0;
}
