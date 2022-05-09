//
// Created by Johan Pedersen on 09/05/2022.
//
#include <iostream>

#ifndef VISIBILITY_2_EXAMPLES_LINKEDLIST_H
#define VISIBILITY_2_EXAMPLES_LINKEDLIST_H

template<class P, class S>
class LinkedList {
private:
    struct Node {
        P data;
        Node *next;
        Node *prev;
        S *top;
        S *bot;

        ~Node() {
            delete top;
            delete bot;
        }
    };

    Node *head;
    Node *tail;


    void deleteList() {

        /* deref head_ref to get the real head */
        Node *current = *head;
        Node *next = nullptr;

        while (current != nullptr) {
            next = current->next;
            delete (current);
            current = next;
        }

        /* deref head_ref to affect the real head back
            in the caller. */
        head = nullptr;
    }

public:
    LinkedList() {
        head = nullptr;
        tail = nullptr;
    }

    ~LinkedList() {
        for(Node* it=head; head; head = head->next){
            delete it;
        }
    }


    void push_front(P d) {
        // Creating new node
        Node *temp;
        temp = new Node();
        temp->data = d;
        temp->prev = nullptr;
        temp->next = head;

        // List is empty
        if (head == nullptr)
            tail = temp;
        else
            head->prev = temp;

        head = temp;
    }

    void insert_before(Node *n, P d) {
        Node *temp;
        temp = new Node();
        temp->data = d;
        temp->next = n;
        temp->prev = n->prev;
        n->prev = temp;

        //if node is to be inserted before first node
        if (n->prev == nullptr)
            head = temp;
    }

    void insert_after(Node *n, P d) {
        Node *temp;
        temp = new Node();
        temp->data = d;
        temp->prev = n;
        temp->next = n->next;
        n->next = temp;

        //if node is to be inserted after last node
        if (n->next == nullptr)
            tail = temp;
    }

    void push_back(P d) {
        // create new node
        Node *temp;
        temp = new Node();
        temp->data = d;
        temp->prev = tail;
        temp->next = nullptr;

        // if list is empty
        if (tail == nullptr)
            head = temp;
        else
            tail->next = temp;
        tail = temp;
    }

    void delete_node(Node *n) {
        // if node to be deleted is first node of list
        if (n->prev == nullptr) {
            head = n->next; //the next node will be front of list
            head->prev = nullptr;
        }
            // if node to be deleted is last node of list
        else if (n->next == nullptr) {
            tail = n->prev;   // the previous node will be last of list
            tail->next = nullptr;
        } else {
            //previous node's next will point to current node's next
            n->prev->next = n->next;
            //next node's prev will point to current node's prev
            n->next->prev = n->prev;
        }
        //delete node
        delete (n);
    }

    void forward_traverse() {
        Node *trav;
        trav = head;
        while (trav != nullptr) {
            std::cout << trav->data << std::endl;
            trav = trav->next;
        }
    }
};

#endif //VISIBILITY_2_EXAMPLES_LINKEDLIST_H
