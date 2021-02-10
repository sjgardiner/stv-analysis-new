#pragma once

// **** Helper code to facilitate setting TTree branch addresses, etc. ****

// A std::unique_ptr with redundant storage of a bare pointer to the managed
// object. This is a hacky workaround for setting TTree branch addresses for
// objects managed by a std::unique_ptr.
template <typename T> class MyPointer : public std::unique_ptr<T> {
  public:

    MyPointer() : std::unique_ptr<T>( new T ) {}

    T*& get_bare_ptr() {
      bare_ptr_ = this->get();
      return bare_ptr_;
    }

  protected:

    T* bare_ptr_ = nullptr;
};

// Helper function template that sets a new address for a pointer to an object
// in an input TTree
template <typename T> void set_object_input_branch_address( TTree& in_tree,
  const std::string& branch_name, T*& address )
{
  in_tree.SetBranchAddress( branch_name.c_str(), &address );
}

// Overloaded version that uses a MyPointer argument instead
template <typename T> void set_object_input_branch_address( TTree& in_tree,
  const std::string& branch_name, MyPointer<T>& u_ptr )
{
  T*& address = u_ptr.get_bare_ptr();
  set_object_input_branch_address( in_tree, branch_name, address );
}

// Helper function that creates a branch (or just sets a new address) for a
// simple variable in an output TTree
void set_output_branch_address( TTree& out_tree, const std::string& branch_name,
  void* address, bool create = false, const std::string& leaf_spec = "" )
{
  if ( create ) {
    if ( leaf_spec != "" ) {
      out_tree.Branch( branch_name.c_str(), address, leaf_spec.c_str() );
    }
    else {
      out_tree.Branch( branch_name.c_str(), address );
    }
  }
  else {
    out_tree.SetBranchAddress( branch_name.c_str(), address );
  }
}

// Helper function template that creates a branch (or just sets a new address)
// for a pointer to an object in an output TTree
template <typename T> void set_object_output_branch_address( TTree& out_tree,
  const std::string& branch_name, T*& address, bool create = false )
{
  if ( create ) out_tree.Branch( branch_name.c_str(), &address );
  else out_tree.SetBranchAddress( branch_name.c_str(), &address );
}

// Overloaded version that uses a MyPointer argument instead
template <typename T> void set_object_output_branch_address( TTree& out_tree,
  const std::string& branch_name, MyPointer<T>& u_ptr, bool create = false )
{
  T*& address = u_ptr.get_bare_ptr();
  set_object_output_branch_address( out_tree, branch_name, address, create );
}
