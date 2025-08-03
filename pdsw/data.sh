mkdir -p data
tar -xzvf dataset_568k_metadata.tar.gz -C data/
make create_ncopy_binary
./create_ncopy_binary ./data/dataset_568k_metadata ./data/dataset_5m_metadata