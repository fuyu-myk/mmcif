use std::sync::Arc;

use crate::document::{Block, Category, Document};

#[derive(Debug, Clone)]
pub struct BlockHandle {
    document: Arc<Document>,
    index: usize,
}

#[derive(Debug, Clone)]
pub struct CategoryHandle {
    block: BlockHandle,
    name: Box<str>,
}

impl BlockHandle {
    pub fn new(document: Arc<Document>, index: usize) -> Option<Self> {
        document.block(index)?;
        Some(BlockHandle { document, index })
    }

    pub fn all(document: &Arc<Document>) -> Vec<BlockHandle> {
        (0..document.len())
            .map(|index| BlockHandle {
                document: Arc::clone(document),
                index,
            })
            .collect()
    }

    pub fn document(&self) -> &Arc<Document> {
        &self.document
    }

    pub fn index(&self) -> usize {
        self.index
    }

    pub fn view(&self) -> Block<'_> {
        self.document
            .block(self.index)
            .expect("block handles are created for existing blocks")
    }

    pub fn category(&self, name: &str) -> Option<CategoryHandle> {
        let found = self.view().category(name)?.name().into();
        Some(CategoryHandle {
            block: self.clone(),
            name: found,
        })
    }
}

impl CategoryHandle {
    pub fn block(&self) -> &BlockHandle {
        &self.block
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn view(&self) -> Category<'_> {
        self.block
            .view()
            .category(&self.name)
            .expect("category handles are created for existing categories")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn handles_keep_the_document_alive() {
        let document =
            Arc::new(Document::parse("data_a\n_x.y 1\ndata_b\nloop_\n_z.w\n1\n2\n").unwrap());
        let blocks = BlockHandle::all(&document);
        drop(document);
        assert_eq!(blocks.len(), 2);
        let category = blocks[1].category("Z").unwrap();
        assert_eq!(category.name(), "z");
        assert_eq!(category.view().row_count(), 2);
        assert_eq!(blocks[0].view().name(), "a");
        assert!(blocks[0].category("z").is_none());
        assert!(BlockHandle::new(Arc::clone(blocks[0].document()), 2).is_none());
    }
}
