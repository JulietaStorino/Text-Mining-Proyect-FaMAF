import tensorflow as tf
from keras.backend.tensorflow_backend import set_session
from neuroner import neuromodel
from shutil import copyfile, copy
import random
import os

#Functions in charge of split the data into an indicate percent
def make_dev_split(data, devPercent):
    index = int(round(devPercent*len(data)))
    shuffled = data[:]
    random.shuffle(shuffled)
    return shuffled[index:], shuffled[:index]

def make_train_valid_files(data_path, valid_percent):
    train_dir = os.path.join(os.getcwd(),"data/train/")
    valid_dir = os.path.join(os.getcwd(),"data/valid/")
    file_names = list(set([os.path.splitext(os.path.join(os.path.abspath(data_path), path))[0] for path in os.listdir(data_path)]))
    train, valid = make_dev_split(file_names, valid_percent)
    for path in train:
        copyfile(path+".ann", train_dir+os.path.basename(path)+".ann")
        copyfile(path+".txt", train_dir+os.path.basename(path)+".txt")    
    for path in valid:
        copyfile(path+".ann", valid_dir+os.path.basename(path)+".ann")
        copyfile(path+".txt", valid_dir+os.path.basename(path)+".txt")

#Function in charge of postprocessing the files resulting from the prediction of the NeuroNER model.
def postprocess(input_file, output_file, clear_ending_points=False, clear_prefixes = False, prefixes = []):
    file = open(output_file, mode='w', encoding='utf8')
    for x in open(input_file).readlines():
        if (len(x)>1):
            token, doc, start, end, e0, e1 = x.split()
            if (clear_ending_points and len(token)>1 and token[-1]=='.'):
                token = token[:-1]
                end = str(int(end)-1)
            if (clear_prefixes and len(token)>1 and token[-1]!='.'):
                for prefix in prefixes:
                    if prefix in token:
                        token = token[len(prefix):]
                        start = str(int(start)+len(prefix))
            file.write(token+" "+doc+" "+start+" "+end+" "+e0+" "+e1+"\n")
        else:
            file.write("\n")
    file.close()

if __name__ == "__main__":
    config = tf.ConfigProto()
    config.gpu_options.allow_growth = True  # dynamically grow the memory used on the GPU
    config.log_device_placement = True  # to log device placement (on which device the operation ran)
    sess = tf.Session(config=config)
    set_session(sess)  # set this TensorFlow session as the default session for Keras

    # Split de datasets
    # Just in case you need to split the dataset into a training set and validation set
    # data_path = './data/'
    # make_train_valid_files(data_path, 0.3)

    # Instantiating and training the NeuroNER model
    #All parameters are taken from the file "parameters.ini". These parameters can be modified when a NeuroNER model is instantiated.
    nn = neuromodel.NeuroNER()
    
    nn.fit()

    # Postprocess task
    run = "old_data/" # the path on which the output data from NeuroNER is allocated
    original_brat = "./data/deploy/"  # the path where we locate the test dataset
    path_to_process = "./output/"+run+"/000_deploy.txt" # the path of the output file on which the entities are predicted by NeuroNER model
    path_to_save = "./output/"+run+"/000_deploy_postprocessed.txt" # the path of the output file on which the post-processed data are saved
    path_to_brat = "./data/deploy/" # the path on which the processed file is going to be saved on brat format
    prefixes = ['nhc/', 'NHC/','nhc:', 'NHC:', 'nhc-', 'NHC-', 'cp:', 'CP:', 'cp-', 'CP-', 'nacimiento:'] # prefixes to be cleaned
    postprocess(path_to_process, path_to_save, clear_ending_points=True, clear_prefixes=True, prefixes=prefixes)

    neuromodel.conll_to_brat.conll_to_brat(path_to_save, path_to_save, original_brat, path_to_brat, overwrite=True)

    #Once the execution has finished, the final output should be placed in the folder indicated by "path_to_brat"

    run = "new_data/" # the path on which the output data from NeuroNER is allocated
    original_brat = "../../data/brat/system/"  # the path where we locate the test dataset
    path_to_process = "./output/"+run+"/000_deploy.txt" # the path of the output file on which the entities are predicted by NeuroNER model
    path_to_save = "./output/"+run+"/000_deploy_postprocessed.txt" # the path of the output file on which the post-processed data are saved
    
    postprocess(path_to_process, path_to_save, clear_ending_points=True, clear_prefixes=True, prefixes=prefixes)

    neuromodel.conll_to_brat.conll_to_brat(path_to_save, path_to_save, original_brat, path_to_brat, overwrite=True)




