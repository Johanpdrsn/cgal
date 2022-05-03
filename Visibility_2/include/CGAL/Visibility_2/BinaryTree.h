template<class T>
class Node {
public:
    T data;
    Node *left;
    Node *right;

    // Val is the key or the value that
    // has to be added to the data part
    explicit Node(T val) {
        data = val;

        // Left and right child for node
        // will be initialized to null
        left = nullptr;
        right = nullptr;
    }

    void printData() {
        for (auto var: data) {
            std::cout << var.vertex(0)->point() << var.vertex(1)->point() << var.vertex(2)->point() << std::endl;
        }
    }
};

template<class T>
class BinaryTree {
public:

    Node<T> *Root() {
        return root;
    }

    void printPreOrder(Node<T> *node) {
        if (node == nullptr)
            return;

        for (auto var: node->data) {
            std::cout << var.vertex(0)->point() << var.vertex(1)->point() << var.vertex(2)->point() << std::endl;
        }
        printPreOrder(node->left);
        printPreOrder(node->right);
    }

    void printInOrder(Node<T> *node) {
        if (node == nullptr)
            return;

        printInOrder(node->left);
        for (auto var: node->data) {
            std::cout << var.vertex(0)->point() << var.vertex(1)->point() << var.vertex(2)->point() << std::endl;
        }
        printInOrder(node->right);
    }

    void printPostOrder(Node<T> *node) {
        if (node == nullptr)
            return;

        printPostOrder(node->left);
        printPostOrder(node->right);
        node->printData();
    }

    void prettyPrint() {
        if (root) {
            printHelper(this->root, "", true);
        }
    }

    BinaryTree(T val) {
        root = new Node<T>(val);
    }

    BinaryTree() {};

private:
    Node<T> *root;

    void printHelper(Node<T> *rootArg, std::string indent, bool last) {
        // print the tree structure on the screen
        if (rootArg != nullptr) {
            std::cout << indent;
            if (last) {
                std::cout << "R----";
                indent += "     ";
            } else {
                std::cout << "L----";
                indent += "|    ";
            }

//            std::cout << rootArg->data.size() << std::endl;
            for (auto var: rootArg->data) {
                std::cout << var.vertex(0)->point() << var.vertex(1)->point() << var.vertex(2)->point() << std::endl;
            }
            printHelper(rootArg->left, indent, false);
            printHelper(rootArg->right, indent, true);
        }
    }
};
