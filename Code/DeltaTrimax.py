
import numpy as np
import sys
import time

class DeltaTrimax():
    def __init__(self,D):
        self.D = D.copy()
        self.D_asli = D.copy()
    
    def hitung_MSR(self, gene, kondisi, waktu, g_add=False, k_add=False, w_add=False):
        
        
        gene_idx = np.expand_dims(np.expand_dims(np.nonzero(gene)[0],axis=0),axis=2)
        #
        waktu_idx = np.expand_dims(np.expand_dims(np.nonzero(waktu)[0],axis=1),axis=1)
        #
        kondisi_idx = np.expand_dims(np.expand_dims(np.nonzero(kondisi)[0],axis=0),axis=0)
        
        
        subarr = self.D[waktu_idx, gene_idx, kondisi_idx]
        self.n_gene = subarr.shape[1]
        self.n_kondisi = subarr.shape[2]
        self.n_waktu = subarr.shape[0]
        
        
        # hitung m_iJK (gene)
        m_iJK = np.nanmean(np.nanmean(subarr,axis=2),axis=0)
        m_iJK = np.expand_dims(np.expand_dims(m_iJK,axis=0),axis=2)
        
        # hitung m_IjK (kondisi)
        m_IjK = np.nanmean(np.nanmean(subarr,axis=0),axis=0)
        m_IjK = np.expand_dims(np.expand_dims(m_IjK,axis=0),axis=1)
        
        # hitung m_IJk (waktu)
        m_IJk = np.nanmean(np.nanmean(subarr,axis=2),axis=1)
        m_IJk = np.expand_dims(np.expand_dims(m_IJk,axis=1),axis=2)
        
        # hitung m_IJK
        m_IJK = np.mean(subarr)
        
        # hitung MSR
        residue = subarr - m_iJK - m_IjK - m_IJk + (2*m_IJK)
        SR = np.square(residue)
        self.MSR = np.mean(SR)
        self.MSR_gene = np.nanmean(np.nanmean(SR,axis=2),axis=0)
        self.MSR_kondisi = np.nanmean(np.nanmean(SR,axis=0),axis=0)
        self.MSR_waktu = np.nanmean(np.nanmean(SR,axis=2),axis=1)
        
        #untuk node addition
        
        if g_add:
            non_gene = np.expand_dims(np.expand_dims(np.nonzero(gene==0)[0],axis=0),axis=2)
            
            
            # hitung m_iJK (untuk gene yang bukan tricluster dari data asli)
            D_b = self.D.copy()
            D_b[waktu_idx,non_gene, kondisi_idx] = self.D_asli[waktu_idx, non_gene, kondisi_idx]

            subarr_b = D_b[waktu_idx, non_gene, kondisi_idx]

            m_iJK_b = np.nanmean(np.nanmean(subarr_b,axis=2),axis=0)
            m_iJK_b = np.expand_dims(np.expand_dims(m_iJK_b,axis=0),axis=2)
            
            r = subarr_b - m_iJK_b - m_IjK - m_IJk + (2*m_IJK)
            sr_b = np.square(r)
            self.MSR_gene_b = np.nanmean(np.nanmean(sr_b,axis=2),axis=0)
            
        
        if k_add:
            non_kondisi = np.expand_dims(np.expand_dims(np.nonzero(kondisi==0)[0],axis=0),axis=0)
            
            # hitung m_IjK (untuk kondisi yg bukan tricluster)
            D_b = self.D.copy()
            D_b[waktu_idx,gene_idx, non_kondisi] = self.D_asli[waktu_idx,gene_idx, non_kondisi]
            subarr_b = D_b[waktu_idx,gene_idx, non_kondisi]

            m_IjK_b = np.nanmean(np.nanmean(subarr_b,axis=0),axis=0)
            m_IjK_b = np.expand_dims(np.expand_dims(m_IjK_b,axis=0),axis=1)
            
            r = subarr_b - m_iJK - m_IjK_b - m_IJk + (2*m_IJK)
            sr_b = np.square(r)
            self.MSR_kondisi_b = np.nanmean(np.nanmean(sr_b,axis=0),axis=0)
        
        if w_add:
            non_waktu = np.expand_dims(np.expand_dims(np.nonzero(waktu==0)[0],axis=1),axis=1)
            
            # hitung m_IJk (untuk waktu yg bukan tricluster)
            D_b = self.D.copy()
            D_b[non_waktu, gene_idx, kondisi_idx] = self.D_asli[non_waktu, gene_idx, kondisi_idx]
            subarr_b = D_b[non_waktu, gene_idx, kondisi_idx]

            m_IJk_b = np.nanmean(np.nanmean(subarr_b,axis=2),axis=1)
            m_IJk_b = np.expand_dims(np.expand_dims(m_IJk_b,axis=1),axis=2)

            r = subarr_b - m_iJK - m_IjK - m_IJk_b + (2*m_IJK)
            sr_b = np.square(r)
            self.MSR_waktu_b = np.nanmean(np.nanmean(sr_b,axis=2),axis=1)
       

    def multiple_node_deletion(self, gene, kondisi, waktu):
        self.hitung_MSR(gene, kondisi, waktu)
        while (self.MSR > self.delta):
            hapus = 0
            
            # hapus gene
            if (self.n_gene > self.gene_cutoff):
                gene_dihapus = self.MSR_gene > (self.MSR * self.lamda)
                nonz_idx = gene.nonzero()[0]
                if gene_dihapus.any():
                    hapus = 1
                gene.put(nonz_idx[gene_dihapus],0)
                # hitung kembali msr
                self.hitung_MSR(gene, kondisi, waktu)
                
            # hapus kondisi
            if (self.n_kondisi > self.kondisi_cutoff):
                kondisi_dihapus = self.MSR_kondisi > (self.MSR * self.lamda)
                nonz_idx = kondisi.nonzero()[0]
                if kondisi_dihapus.any():
                    hapus = 1
                kondisi.put(nonz_idx[kondisi_dihapus],0)
                # hitung kembali msr
                self.hitung_MSR(gene, kondisi, waktu)
                
            # hapus waktu
            if (self.n_waktu > self.waktu_cutoff):
                waktu_dihapus = self.MSR_waktu > (self.MSR * self.lamda)
                nonz_idx = waktu.nonzero()[0]
                if waktu_dihapus.any():
                    hapus = 1
                waktu.put(nonz_idx[waktu_dihapus],0)
                # hitung kembali msr
                self.hitung_MSR(gene, kondisi, waktu)
            
            # menghentikan iterasi
            if not hapus:
                break
            sys.stdout.write('\r multiple deletion {} - single deletion {} - node additon {}'.format(
                self.muldel, self.singdel, self.nodeadd
            ))
            sys.stdout.flush()

            self.muldel += 1

            
        return gene, kondisi, waktu
  
    def single_node_deletion(self, gene, kondisi, waktu):
        self.hitung_MSR(gene, kondisi, waktu)
        while (self.MSR > self.delta):
            
            gene_max = np.argmax(self.MSR_gene)
            kondisi_max = np.argmax(self.MSR_kondisi)
            waktu_max = np.argmax(self.MSR_waktu)
            
            if (self.MSR_gene[gene_max] > self.MSR_kondisi[kondisi_max]):
                if (self.MSR_gene[gene_max] > self.MSR_waktu[waktu_max]):
                    # Hapus Gene
                    nonz_idx = gene.nonzero()[0]
                    gene.put(nonz_idx[gene_max],0)
                else:
                    # Hapus waktu
                    nonz_idx = waktu.nonzero()[0]
                    waktu.put(nonz_idx[waktu_max],0)
                    #print("menghapus waktu ke-",waktu_max)
            else:
                if (self.MSR_kondisi[kondisi_max] > self.MSR_waktu[waktu_max]):
                    # Hapus Kondisi
                    nonz_idx = kondisi.nonzero()[0]
                    kondisi.put(nonz_idx[kondisi_max],0)
                    #print("menghapus kondisi ke-",kondisi_max)
                else:
                    # Hapus Waktu
                    nonz_idx = waktu.nonzero()[0]
                    waktu.put(nonz_idx[waktu_max],0)
                    #print("menghapus waktu ke-",waktu_max)
            
            sys.stdout.write('\r multiple deletion {} - single deletion {} - node additon {}'.format(
                self.muldel, self.singdel, self.nodeadd
            ))
            sys.stdout.flush()
            
            #i+=1
            self.singdel += 1
            self.hitung_MSR(gene, kondisi, waktu)


        #print("MSR akhir single node deletion :", self.MSR)
        return gene, kondisi, waktu
                    
    def node_addition(self, gene, kondisi, waktu):
        self.hitung_MSR(gene, kondisi, waktu)
        while True:
            #print("node addition ke-",i)
            
            self.hitung_MSR(gene, kondisi, waktu)
            
            #penambaha gen
            self.hitung_MSR(gene, kondisi, waktu, g_add=True)
            no_gene_idx = np.nonzero(gene==0)[0]
            gene_to_add = self.MSR_gene_b < self.MSR
            if gene_to_add.any():
              gene.put(no_gene_idx[gene_to_add],1)

              g_idx = np.expand_dims(np.expand_dims(no_gene_idx[gene_to_add],axis=0),axis=2)
              k_idx = np.expand_dims(np.expand_dims(np.nonzero(kondisi)[0],axis=0),axis=0)
              w_idx = np.expand_dims(np.expand_dims(np.nonzero(waktu)[0],axis=1),axis=2)

              self.D[w_idx, g_idx, k_idx] = self.D_asli[w_idx, g_idx, k_idx]
             
            self.hitung_MSR(gene, kondisi, waktu)
             
            #penambahan kondisi
            self.hitung_MSR(gene, kondisi, waktu, k_add =True)
            no_kondisi = np.nonzero(kondisi==0)[0]
            kondisi_to_add = self.MSR_kondisi_b < self.MSR
            if kondisi_to_add.any():
              kondisi.put(no_kondisi[kondisi_to_add],1)

              g_idx = np.expand_dims(np.expand_dims(np.nonzero(gene)[0],axis=0),axis=2)
              k_idx = np.expand_dims(np.expand_dims(no_kondisi[kondisi_to_add],axis=0),axis=0)
              w_idx = np.expand_dims(np.expand_dims(np.nonzero(waktu)[0],axis=1),axis=2)

              self.D[w_idx, g_idx, k_idx] = self.D_asli[w_idx, g_idx, k_idx]
             
            self.hitung_MSR(gene, kondisi, waktu)
             
            #penambaha waktu
            self.hitung_MSR(gene, kondisi, waktu, w_add=True)
            no_waktu = np.nonzero(waktu==0)[0]
            waktu_to_add = self.MSR_waktu_b < self.MSR
            if waktu_to_add.any():
              waktu.put(no_waktu[waktu_to_add],1)

              g_idx = np.expand_dims(np.expand_dims(np.nonzero(gene)[0],axis=0),axis=2)
              k_idx = np.expand_dims(np.expand_dims(np.nonzero(kondisi)[0],axis=0),axis=0)
              w_idx = np.expand_dims(np.expand_dims(no_waktu[waktu_to_add],axis=1),axis=2)

              self.D[w_idx, g_idx, k_idx] = self.D_asli[w_idx, g_idx, k_idx]
             
            self.hitung_MSR(gene, kondisi, waktu)

            if not gene_to_add.any() and not kondisi_to_add.any() and not waktu_to_add.any():
                break
            
            sys.stdout.write('\r multiple deletion {} - single deletion {} - node additon {}'.format(
                self.muldel, self.singdel, self.nodeadd
            ))
            sys.stdout.flush()

            #i+=1
            self.nodeadd +=1
        #print("MSR akhir node addition :",self.MSR)
        return gene, kondisi, waktu
    
    
    def mask(self, gene, kondisi, waktu, minval, maxval):
        g = np.expand_dims(np.expand_dims(gene.nonzero()[0],axis=0),axis=2)
        k = np.expand_dims(np.expand_dims(kondisi.nonzero()[0],axis=0),axis=0)
        w = np.expand_dims(np.expand_dims(waktu.nonzero()[0],axis=1),axis=2)
        
        shape = np.count_nonzero(waktu), np.count_nonzero(gene), np.count_nonzero(kondisi)
        mask_val = np.random.uniform(self.minval, self.maxval, shape)
        self.D[w,g,k] = mask_val
    
    def TQI(self,msr,gene,kondisi,waktu):
      	m = self.D_asli[waktu][:,gene][:,:,kondisi]
      	x,y,z = m.shape
      	v = x*y*z
      	tqi = msr/v
      	return tqi
            

    
    def fit(self, delta, lamda, n_triclusters=0):

        awal = time.time()
        
        n_waktu, n_gene, n_kondisi = self.D.shape
        
        self.delta = delta
        self.lamda = lamda
        
        # treshold untuk multiple node deletion
        self.gene_cutoff, self.kondisi_cutoff, self.waktu_cutoff = 50,50,50
        #nilai untuk masking
        self.minval = np.min(self.D)
        self.maxval = np.max(self.D)
        
        hasil_gen = []
        hasil_kondisi = []
        hasil_waktu = []
        
        msr = []
        TQI =[]

        i = 1
        
        while True:
            self.muldel = 0
            self.singdel =0
            self.nodeadd = 0

            print("Tricluster ",i)
            waktu = np.ones(n_waktu, dtype=np.bool)
            gene = np.ones(n_gene, dtype=np.bool)
            kondisi = np.ones(n_kondisi, dtype=np.bool)
            
            # Multiple Node Deletion
            gene, kondisi, waktu = self.multiple_node_deletion(gene, kondisi, waktu)
            
            # Single Node Deletion
            gene, kondisi, waktu = self.single_node_deletion(gene, kondisi, waktu)
            
            # Node Addition
            gene, kondisi, waktu = self.node_addition(gene, kondisi, waktu)
            
            # TQI
            tqi = self.TQI(self.MSR, gene, kondisi, waktu )
            

            if (gene.sum()==1) or (kondisi.sum()==1) or (waktu.sum()==1):
                akhir = ((time.time()-awal)/60)
                print("\n Waktu Komputasi : {} menit".format(akhir))
                break
            
            print("\n------------- MSR: ",self.MSR)
            print("------------- TQI: ",tqi)

            hasil_gen.append(gene)
            hasil_kondisi.append(kondisi)
            hasil_waktu.append(waktu)

            msr.append(self.MSR)
            TQI.append(tqi)


            if (n_triclusters == i):
                break
    
    
            # mask
            self.mask(gene, kondisi, waktu, self.minval, self.maxval)
            
            i+=1
            
        return hasil_gen, hasil_kondisi, hasil_waktu, msr, tqi