"""Create the evidence-first schema.

Revision ID: 20260817_0001
Revises: None
"""

from alembic import op

from src.workbench.models import Base

revision = "20260817_0001"
down_revision = None
branch_labels = None
depends_on = None


def upgrade() -> None:
    # The initial migration is intentionally derived from the reviewed ORM
    # metadata so SQLCipher and test SQLite receive identical constraints.
    Base.metadata.create_all(bind=op.get_bind())


def downgrade() -> None:
    Base.metadata.drop_all(bind=op.get_bind())
