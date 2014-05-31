#ifndef __COUNT_PTR_HPP__
#define __COUNT_PTR_HPP__


#include <cstddef>            // std::::std::size_t
#include <algorithm>          // std::swap

#if 0
#include <boost/shared_ptr.hpp>
#include <tr1/shared_ptr.hpp> 
#endif

namespace mallie{
	
	template<class T>
	class direct_count_ptr{
	public:
		class counter{
		public:
			explicit counter():ctr(1){}	
			~counter()throw(){}
			
			void up() throw(){ctr++;}
			void down() throw(){ctr--;}
			
			int  count()const{return ctr;}	
			bool is_only()const throw(){return (ctr < 2);}
			
		protected:
			int ctr;		//counter	
		};
	public:
		typedef direct_count_ptr<T> this_type;
		typedef T element_type;
	public:
		explicit direct_count_ptr(T * p = 0):ptr(p),ctr(new counter()){}
		
		direct_count_ptr(const this_type & rhs):ptr(rhs.ptr),ctr(rhs.ctr){
			ctr->up();	
		}
		
		~direct_count_ptr() throw() {		
			if(ctr->is_only()){
				delete ptr;
				delete ctr;
			}else{
				ctr->down();	
			}
		}
		
		this_type& operator=(const this_type & rhs){
			if(&rhs == this)return *this;
			
			if(ctr->is_only()){
				delete ptr;
				delete ctr;
			}else{
				ctr->down();	
			}
			ptr = rhs.ptr;
			ctr = rhs.ctr;
			ctr->up();
			
			return *this;
		}
		
		void reset(T* p = 0){
			if(p != ptr){
				if(ctr->is_only()){
					delete ptr;
					delete ctr;
				}else{
					ctr->down();	
				}
				ptr = p;
				ctr = new counter(); 
			}
		}
		
		bool is_only()const throw(){ return ctr->is_only();	}	
		T* operator-> () const throw(){	return ptr; }
		T* get() const throw(){ return ptr; }	
		T& operator* () const throw(){ return *ptr;}
		
		T* get_pointer() const throw(){ return ptr; }//boost::mem_fn	
		
		void swap( this_type & other ) throw() {
			std::swap(ptr,other.ptr);
			std::swap(ctr,other.ctr);
		}
		
	private:
		counter* ctr;
		T* ptr;
	};
	

	template<class T>
	class count_ptr{
	public:
		class counter_with_pointer{
		public:
			explicit counter_with_pointer(T* p = 0):ctr(1),ptr(p){}	
			~counter_with_pointer(){if(ptr != 0){delete ptr;}}
			
			T* get() const throw(){return ptr;}	
			T& operator* () const throw(){return *ptr;}	
			T* operator-> () const throw(){return ptr;}

			void up() throw(){ctr++;}
			void down() throw(){ctr--;}

			int  count()const{return ctr;}	
			bool is_only()const throw(){return (ctr < 2);}
			
		private:
			int ctr;
			T * ptr;		//inside pointer			
		};
		typedef counter_with_pointer rc_type;
	public:
		typedef count_ptr<T> this_type;
		typedef T element_type;
	public:
		explicit count_ptr(T * p = 0):ptr(new rc_type(p)){}
		
		count_ptr(const this_type & rhs):ptr(rhs.ptr){
			ptr->up();	
		}
		
		~count_ptr() throw() {		
			if(ptr->is_only()){
				delete ptr;
			}else{
				ptr->down();	
			}
		}
		
		this_type& operator=(const this_type & rhs){
			if(&rhs == this)return *this;

			if(ptr->is_only()){
				delete ptr;
			}else{
				ptr->down();	
			}
			ptr = rhs.ptr;
			ptr->up();
			
			return *this;
		}
		
		void reset(T* p = 0){
			if(p != ptr->get()){
				if(ptr->is_only()){
					delete ptr;
				}else{
					ptr->down();	
				}
				ptr = new rc_type(p);
			}
		}
		
		bool is_only()const throw(){ return ptr->is_only();	}	
		T* operator-> () const throw(){	return ptr->operator->(); }
		T* get() const throw(){ return ptr->get(); }	
		T& operator* () const throw(){ return *(ptr->get());}

    int  count()const{return ptr->count();}
    long use_count()const{return ptr->count();}//shared_ptr compatible
		
		void swap( this_type & other ) throw() {
			std::swap(ptr,other.ptr);		
		}
		
	private:
		rc_type* ptr;
	};	

	template<class T>
	class auto_count_ptr:public count_ptr<T>{
	public:
		auto_count_ptr(T* ptr):	count_ptr<T>(ptr){}
		auto_count_ptr(const count_ptr<T>& ptr):count_ptr<T>(ptr){}
		auto_count_ptr(const auto_count_ptr<T>& ptr):count_ptr<T>(ptr){}
	};

}

namespace std{
	template<class T>
	inline void swap(mallie::count_ptr<T> & lhs,  mallie::count_ptr<T> & rhs){
		lhs.swap(rhs);
	}
	
	template<class T>
	inline void swap(mallie::direct_count_ptr<T> & lhs, mallie::direct_count_ptr<T> & rhs){
		lhs.swap(rhs);
	}
}

#endif

				   
