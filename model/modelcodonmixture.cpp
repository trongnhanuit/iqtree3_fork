//
//  modelcodonmixture.cpp
//  iqtree
//
//  Created by Minh Bui on 4/3/2025.
//

#include "modelcodonmixture.h"


ModelCodonMixture::ModelCodonMixture(string orig_model_name, string model_name,
    ModelsBlock *models_block, StateFreqType freq, string freq_params,
    PhyloTree *tree, bool optimize_weights)
: ModelMarkov(tree), ModelMixture(tree)
{
    if (tree->aln->seq_type != SEQ_CODON)
        outError("Can't apply codon mixture model as sequence type is not codon");
    auto cmix_pos = orig_model_name.find("+CMIX");
    ASSERT(cmix_pos != string::npos);
    auto end_pos = orig_model_name.find_first_of("+*{", cmix_pos+1);
    string cmix_type;
    if (end_pos == string::npos)
        cmix_type = orig_model_name.substr(cmix_pos+5);
    else
        cmix_type = orig_model_name.substr(cmix_pos+5, end_pos-cmix_pos-5);
    
    // read the input parameters for +CMIXi{...}
    StrVector vec;
    if (end_pos != string::npos && orig_model_name[end_pos] == '{') {
        auto close_br_pos = orig_model_name.find_first_of("}", end_pos+1);
        if (close_br_pos != string::npos) {
            string cmix_param_list = orig_model_name.substr(end_pos+1,close_br_pos-end_pos-1);
            convert_string_vec(cmix_param_list.c_str(), vec);
        }
    }
    
    string model_list = "";
    /* setting fix_kappa for class 2 and 3 */
    if (vec.size() == 0) {
        if (cmix_type == "1a") {
            // M1a neural model with 2 classes, omega2 = 1.0
            model_list = model_name + "{<0.999}," + model_name + "{1.0,1.0}";
        } else if (cmix_type == "2a") {
            // M2a selection model with 3 classes
            model_list = model_name + "{<0.999}," + model_name + "{1.0,1.0}," + model_name + "{>1.001,1.0}";
        } else if (cmix_type == "3") {
            // M3 model with 3 classes with no constraint
            model_list = model_name + "{>0.001}," + model_name + "{>0.001,1.0}," + model_name + "{>0.001,1.0}";
        } else {
            outError("Unknown codon mixture " + orig_model_name.substr(cmix_pos));
        }
    } else {
        // user inputs the parameter values for CMIX model
        if (cmix_type == "1a" && vec.size() != 2 && vec.size() != 4)
            outError("Error! There should be 2 (or 4) parameters inside CMIX1a{} stating the omega value (and the weight) of each class");
        else if (cmix_type == "2a" && vec.size() != 3 && vec.size() != 6)
            outError("Error! There should be 3 (or 6) parameters inside CMIX2a{} stating the omega value (and the weight) of each class.");
        else if (cmix_type == "3" && vec.size() != 3 && vec.size() != 6)
            outError("Error! There should be 3 (or 6) parameters inside CMIX3{} stating the omega value (and the weight) of each class.");
        if (vec.size() >= 4) {
            // with both omega and weight
            for (int i = 0; i < vec.size(); i+=2) {
                if (i > 0) {
                    model_list += ",";
                }
                model_list += model_name + "{" + vec[i] + ",1.0}:1.0:" + vec[i+1];
            }
        } else {
            // only omega
            for (int i = 0; i < vec.size(); i++) {
                if (i > 0) {
                    model_list += ",";
                }
                model_list += model_name + "{" + vec[i] + ",1.0}";
            }
        }
    }
    initMixture(orig_model_name, model_name, model_list, models_block,
                freq, freq_params, tree, optimize_weights);
    
    // show the initial parameters
    cout << "Initial parameters in the Codon Mixture:" << endl;
    writeInfo(cout);
    cout << endl;
}

ModelCodonMixture::~ModelCodonMixture()
{
    
}

bool ModelCodonMixture::getVariables(double *variables) {
    bool changed = ModelMixture::getVariables(variables);
    auto kappa = ((ModelCodon*)at(0))->kappa;
    auto kappa2 = ((ModelCodon*)at(0))->kappa2;
    for (int i = 1; i < size(); i++) {
        ModelCodon *model = (ModelCodon*)at(i);
        model->kappa = kappa;
        model->kappa2 = kappa2;
    }
    return changed;
}


void ModelCodonMixture::setVariables(double *variables) {
    ModelMixture::setVariables(variables);
}
